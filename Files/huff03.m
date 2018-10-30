function varargout = huff03(arg1, arg2, arg3)
% Huffman encoder/decoder, one or two vectors of non-negative
% integers are Huffman coded.
% Note that number of arguments decide whether it is encoding og decoding
% If it is more than one input, encoding is done
%
% [y, br, bre] = huff03(x1, x2, Speed);       % encoding
% [y, br] = huff03(x1, x2);                   % encoding
% y = huff03(x1,0);                           % encoding
% x1 = huff03(y);                             % decoding
% [x1, x2] = huff03(y);                       % decoding
% ------------------------------------------------------------------
% Arguments:
%  x1, x2   a column vector of non-negative integers representing the
%           symbol sequence.
%           Best results are achieved if x1, x2 are mostly small integers
%  y        a column vector of non-negative integers (bytes) representing 
%           the code, 0 <= y(i) <= 255. Note that the first part of y will
%           be the Huffman-table/tree information.
%  Speed    Set this to 1 to cheat during encoding.
%           y will then be a sequence of zeros only, but it will be
%           of correct length and the other output arguments will 
%           be correct. If last argumet is scalar, then it is assumed
%           to be 'Speed'.
% Optional output arguments
%  br       Actual bit rate, br(1) for x1, br(2) for x2, br(3) overall
%  bre      0th order entropy (minimum bit rate) of the symbols. 
%           bre(1) is for x1, bre(2) is for x2, bre(3) is overall
% ------------------------------------------------------------------

% SOME NOTES ON THE FUNCTION
% huff03 depens on other functions for Huffman code, and the functions in this file
%  hufflen  - find length of codewords
%  huffcode - find huffman codewords
%  hufftree - find huffman tree

%----------------------------------------------------------------------
% Copyright (c) 1999.  Karl Skretting.  All rights reserved.
% Hogskolen in Stavanger (Stavanger University), Signal Processing Group
% Mail:  karl.skretting@tn.his.no   Homepage:  http://www.ux.his.no/~karlsk/
% 
% HISTORY:
% Ver. 1.0  17.01.99  KS: Function made as part of Signal Compression Project 98
% Ver. 1.1  19.01.99  KS: added Level and TryAll to get better compression
%                     at the cost of more time and complexity
%----------------------------------------------------------------------

global y Byte BitPos Speed Level TryAll

Level=0;     % start level
TryAll=1;    % 1 - best compression   0 - fastest compression

% check input and output arguments
if (nargin < 1); 
   error('huff03: function must have input arguments, see help.'); 
end
if (nargout < 1); 
   error('huff03: function must have output arguments, see help.'); 
end
if (nargin == 1)
   Encode=0;Decode=1;
else
   Encode=1;Decode=0;
end;   

% assigne values to arguments
if Decode
   y=arg1(:);
end
if Encode
   x1=arg1(:);
   if ((nargin==2) & (length(arg2(:))==1))
      Speed=arg2;x2=[];
   elseif ((nargin==2) & (length(arg2(:))>1))
      x2=arg2(:);
      Speed=0;
   elseif ((nargin==3) & (length(arg3(:))==1))
      x2=arg2(:);
      Speed=arg3;
   else
      error('huff03: wrong input arguments, see help.'); 
   end
end
clear arg1 arg2 arg3
   
if Encode
   L1=length(x1);L2=length(x2);
   br=[0,0,0];bre=[0,0,0];
   % initalize the global variables
   y=zeros(2*(L1+L2),1);
   Byte=0;BitPos=1;  % ready to write into first position
   % start encoding, first bit indicate one or two sequences
   if (L2 == 0) 
      if Speed
         BitPos=BitPos-1;
         if (~BitPos); Byte=Byte+1; BitPos=8; end; 
      else
         PutBit(0); 
      end;
      [bits1, ent1]=EncodeVector(x1);
      ent2=0;bits2=0;
   else
      if Speed
         BitPos=BitPos-1;
         if (~BitPos); Byte=Byte+1; BitPos=8; end; 
      else
         PutBit(1); 
      end;
      [bits1, ent1]=EncodeVector(x1);
      [bits2, ent2]=EncodeVector(x2);
   end
   y=y(1:Byte);   
   varargout(1) = {y};
   if (nargout >= 2)
      if (L1>0); br(1)=bits1/L1; end;
      if (L2>0); br(2)=bits2/L2; end;
      br(3)=(1+bits1+bits2)/(L1+L2);
      varargout(2) = {br};
   end
   if (nargout >= 3)
      bre(1)=ent1;
      bre(2)=ent2;
      bre(3)=(L1*ent1+L2*ent2)/(L1+L2);
      varargout(3) = {bre};
   end
end
if Decode
   % initalize the global variables, y is set earlier
   Byte=0;BitPos=1;  % ready to read from first position
   % start decoding, is it one or two sequences
   if GetBit         
      x1=DecodeVector;
      x2=DecodeVector;
   else
      x1=DecodeVector;
      x2=[];
   end
   if (nargout == 1)
      if (length(x2)>0)
         x1=[x1;x2];
         warning('huff03: two vectors are concatenated.');
      end
      varargout(1) = {x1};
   end
   if (nargout == 2)
      varargout(1) = {x1};
      varargout(2) = {x2};
   end
end

return     % end of main function, huff03

% the EncodeVector and DecodeVector functions are the ones
% where actual coding is going on.
function [bits, ent] = EncodeVector(x, bits, L, HL, Method, Maxx, Meanx)
global y Byte BitPos Speed Level TryAll
Level = Level + 1;
if (nargin==1)            % this is the same as (Level==1)
   L=length(x);
   Maxx=max(x);
   Meanx=mean(x);
   % find the histogram, hist function is slow if many bins
   if (Maxx < 72)
      Hi=hist(x,0:Maxx);
   else
      Hi=zeros(Maxx+1,1);
      for l=1:L; Hi(x(l)+1)=Hi(x(l)+1)+1; end;
   end
   Hinz=nonzeros(Hi);
   ent=log2(L)-sum(Hinz.*log2(Hinz))/L;  % find entropy
   HL=hufflen(Hi);
   if ((Maxx>250) | TryAll)
      HLlen=HuffTabLen(HL,0);
   else
      HLlen=HuffTabLen(HL,1);
   end   
   I=find(HLlen==min(HLlen));Method=I(1);
   % find number of bits to use
   bits=6;
   bits=bits+HLlen(Method);
   bits=bits+sum(HL.*Hi);
   if ((Maxx<250) & TryAll)
      % test if 'non-optimal' Huffman Code would give fewer bits 
      % (due to shorte table), this makes the function a little bit slower
      HL1=hufflen(Hi+1);
      HLlen1=HuffTabLen(HL1,0);
      I1=find(HLlen1==min(HLlen1));Method1=I1(1);
      % find number of bits to use
      bits1=6;
      bits1=bits1+HLlen1(Method1);
      bits1=bits1+sum(HL1.*Hi);
      if (bits1 < bits)               % Use these codeword lengths instead
         bits=bits1;Method=Method1;
         HL=HL1;HLlen=HLlen1;
      end
   end
   if (L>=16); bits=bits+4; end;
   if (L>=272); bits=bits+4; end;
   if (L>=4368); bits=bits+4; end;
else                % arguments are given, do not need to be calculated
   ent=0;
end
%
% Here we have: x, bits, L, HL, Method, Maxx, Meanx, ent
if ((Level < 8) & TryAll & (L>10))           % maximum level is set to 8
   xm=median(x);
   % xm=mean(x);
   x1=zeros(L,1);x2=zeros(L,1);
   x2(1)=x(1);i1=0;i2=1;
   for i=2:L
      if (x(i-1) <= xm) 
         i1=i1+1; x1(i1)=x(i);
      else
         i2=i2+1; x2(i2)=x(i);
      end
   end
   x1=x1(1:i1);x2=x2(1:i2);
   % find bits1 and bits2 for x1 and x2
   L1=length(x1);L2=length(x2);
   Maxx1=max(x1);Maxx2=max(x2);
   Meanx1=mean(x1);Meanx2=mean(x2);
   if (Maxx1 < 72)
      Hi1=hist(x1,0:Maxx1);
   else
      Hi1=zeros(Maxx1+1,1);
      for l=1:L1; Hi1(x1(l)+1)=Hi1(x1(l)+1)+1; end;
   end
   if (Maxx2 < 72)
      Hi2=hist(x2,0:Maxx2);
   else
      Hi2=zeros(Maxx2+1,1);
      for l=1:L2; Hi2(x2(l)+1)=Hi2(x2(l)+1)+1; end;
   end
   HL1=hufflen(Hi1);HL2=hufflen(Hi2);
   if ((Maxx1>250) | TryAll); HLlen1=HuffTabLen(HL1,0); else; HLlen1=HuffTabLen(HL1,1); end;
   if ((Maxx2>250) | TryAll); HLlen2=HuffTabLen(HL2,0); else; HLlen2=HuffTabLen(HL2,1); end;
   I=find(HLlen1==min(HLlen1));Method1=I(1);
   I=find(HLlen2==min(HLlen2));Method2=I(1);
   bits1=6+HLlen1(Method1)+sum(HL1.*Hi1);
   bits2=6+HLlen2(Method2)+sum(HL2.*Hi2);
   if ((Maxx1<250) & TryAll)
      HL3=hufflen(Hi1+1);
      HLlen3=HuffTabLen(HL3,0);
      I=find(HLlen3==min(HLlen3));Method3=I(1);
      bits3=6+HLlen3(Method3)+sum(HL3.*Hi1);
      if (bits3 < bits1); bits1=bits3;Method1=Method3;HL1=HL3;HLlen1=HLlen3; end;
   end
   if ((Maxx2<250) & TryAll)
      HL3=hufflen(Hi2+1);
      HLlen3=HuffTabLen(HL3,0);
      I=find(HLlen3==min(HLlen3));Method3=I(1);
      bits3=6+HLlen3(Method3)+sum(HL3.*Hi2);
      if (bits3 < bits2); bits2=bits3;Method2=Method3;HL2=HL3;HLlen2=HLlen3; end;
   end
   if (L1>=16); bits1=bits1+4; end;
   if (L1>=272); bits1=bits1+4; end;
   if (L1>=4368); bits1=bits1+4; end;
   if (L2>=16); bits2=bits2+4; end;
   if (L2>=272); bits2=bits2+4; end;
   if (L2>=4368); bits2=bits2+4; end;
else
   bits1=bits;bits2=bits;
end
% Here we may have: x1, bits1, L1, HL1, Method1, Maxx1, Meanx1
% and               x2, bits2, L2, HL2, Method2, Maxx2, Meanx2
% but at least we have bits1 and bits2
if ((bits1+bits2) < bits)
   if Speed
      BitPos=BitPos-1;
      if (~BitPos); Byte=Byte+1; BitPos=8; end; 
   else
      PutBit(1); 
   end;
   [bits1, temp] = EncodeVector(x1, bits1, L1, HL1, Method1, Maxx1, Meanx1);
   [bits1, temp] = EncodeVector(x2, bits2, L2, HL2, Method2, Maxx2, Meanx2);
   bits=bits1+bits2+1;
else
   bits=bits+1; 
   % disp(['huff03-EncodeVector: Level=',int2str(Level),'  ',int2str(L),...
   %       ' sybols stored in ',int2str(bits),' bits.']);
   if Speed
      % advance Byte and BitPos without writing to y
      Byte=Byte+floor(bits/8);
      BitPos=BitPos-mod(bits,8);
      if (BitPos<=0); BitPos=BitPos+8; Byte=Byte+1; end;
   else
      % put the bits into y
      StartPos=Byte*8-BitPos;
      PutBit(0);
      PutVLIC(L);        % max for L is 69903
      PutHuffTab(HL,Method);
      HK=huffcode(HL);
      for i=1:L;
         n=x(i)+1;    % symbol number (value 0 is first symbol, symbol 1)
         for k=1:HL(n)
            PutBit(HK(n,k));
         end
      end
      % check if one has used as many bits as calculated
      BitsUsed=Byte*8-BitPos-StartPos;
      if (BitsUsed~=bits)
         disp(['BitsUsed=',int2str(BitsUsed),'  bits=',int2str(bits)]);
         error('huff03-EncodeVector: Logical error, (BitsUsed~=bits).'); 
      end
   end
end
Level = Level - 1;
return    % end of EncodeVector

function x = DecodeVector
global y Byte BitPos
if GetBit
   x1=DecodeVector;
   x2=DecodeVector;
   L=length(x1)+length(x2);
   xm=median([x1;x2]);
   % xm=mean([x1;x2]);
   x=zeros(L,1);
   x(1)=x2(1);
   i1=0;i2=1;
   for i=2:L
      if (x(i-1) <= xm) 
         i1=i1+1; x(i)=x1(i1);
      else
         i2=i2+1; x(i)=x2(i2);
      end
   end
else
   L=GetVLIC;
   x=zeros(L,1);
   HL=GetHuffTab;
   Htree=hufftree(HL);
   root=1;pos=root;     
   l=0;  % number of symbols decoded so far
   while l<L
      if GetBit
         pos=Htree(pos,3);
      else
         pos=Htree(pos,2);
      end
      if Htree(pos,1)           % we have arrived at a leaf
         l=l+1;
         x(l)=Htree(pos,2)-1;   % value is one less than symbol number
         pos=root;              % start at root again
      end
   end
end
return    % end of DecodeVector

% Functions to write and read the Huffman Table Information
% Several possible schemes may be used to store the 
% Huffman Codeword Lengths, normally the best one is
% chosen. For decoding we do not need to know the method
% since it is always include in the two first bytes
% HuffTabLen returns how many bit needed to store Huffman-table/tree 
% Method is from 1 to 4, if Method is 0, all metheods are
% done and HLlen is returned as a vector
%  Method 1: (max codeword length is 15)
%         2 bit for method, '00'
%         6-18 bit is VLIC for number of symbols
%         4 bit to give length of first symbol
%        Then for each of the next symbols a code to tell its length
%         '0'             - same length as previous symbol
%         '10'            - increase length by 1
%         '1100'          - reduce length by 1
%         '1101'          - increase length by 2
%         '111xxxx'       - set symbol length to xxxx
%  Method 2: (max codeword length is 15)
%         2 bit for method, '01'
%         6-18 bit is VLIC for number of symbols
%         4 bit to give length of each symbol
%  Method 3: (max codeword length is 15)
%         2 bit for method, '10'
%         6-18 bit is VLIC for number of symbols
%         4 bit to give length of first symbol
%        Then for next symbols
%         '00'+VLIC       - VLIC symbols with same length as previous symbol
%         '01'            - increase length by 1
%         '10'            - reduce length by 1
%         '110'           - increase length by 2
%         '111xxxx'       - set symbol length to xxxx
%  Method 4: (max codeword length is 31)
%         2 bit for method, '11'
%         6-18 bit is VLIC for number of symbols
%         5 bit to give length of first symbol
%        Then for next symbols
%         '0'+VLIC        - VLIC symbols with same length as previous symbol
%         '10'            - increase length by 1
%         '11xxxxx'       - set symbol length to xxxxx
%
function HLlen = HuffTabLen(HL,Method)
L=length(HL);
StartBits=8;
if (L>=16); StartBits=StartBits+4; end;
if (L>=272); StartBits=StartBits+4; end;
if (L>=4368); StartBits=StartBits+4; end;
if (max(HL) >= 16); Method=4; end;
if (Method==0); HLlen=[1,1,1,1]*999999; else; HLlen=999999; end;
HLi=0;
if ((Method==0) | (Method==1))
   b=StartBits+4;    % include first length
   for l=2:L
      if (HL(l)==HL(l-1)); b=b+1;
      elseif (HL(l)==(HL(l-1)+1)); b=b+2;
      elseif (HL(l)==(HL(l-1)-1)); b=b+4;
      elseif (HL(l)==(HL(l-1)+2)); b=b+4;
      else; b=b+7;
      end
   end
   HLi=HLi+1;
   HLlen(HLi)=b;   
end
if ((Method==0) | (Method==2))
   HLi=HLi+1;
   HLlen(HLi)=StartBits+4*L;   
end
if ((Method==0) | (Method==3))
   b=StartBits+4;    % include first length
   r=0;
   for l=2:L
      if (HL(l)==HL(l-1)); r=r+1;
      else
         if r>0
            b=b+8;
            if (r>=16); b=b+4; end;
            if (r>=272); b=b+4; end;
            if (r>=4368); b=b+4; end;
            r=0;
         end
         if (HL(l)==(HL(l-1)+1)); b=b+2; 
         elseif (HL(l)==(HL(l-1)-1)); b=b+2;
         elseif (HL(l)==(HL(l-1)+2)); b=b+3;
         else b=b+7;
         end
      end
   end
   if r>0
      b=b+8;
      if (r>=16); b=b+4; end;
      if (r>=272); b=b+4; end;
      if (r>=4368); b=b+4; end;
      r=0;
   end
   HLi=HLi+1;
   HLlen(HLi)=b;   
end
if ((Method==0) | (Method==4))
   b=StartBits+5;    % include first length
   r=0;
   for l=2:L
      if (HL(l)==HL(l-1)); r=r+1;
      else
         if r>0
            b=b+7;
            if (r>=16); b=b+4; end;
            if (r>=272); b=b+4; end;
            if (r>=4368); b=b+4; end;
            r=0;
         end
         if (HL(l)==(HL(l-1)+1)); b=b+2; 
         else b=b+7;
         end
      end
   end
   if r>0
      b=b+7;
      if (r>=16); b=b+4; end;
      if (r>=272); b=b+4; end;
      if (r>=4368); b=b+4; end;
      r=0;
   end
   HLi=HLi+1;
   HLlen(HLi)=b;   
end
return;  % end of HuffTabLen

function PutHuffTab(HL,Method)
global y Byte BitPos
L=length(HL);
if (max(HL) >= 16); Method=4; end;
if (Method==1)
   PutBit(0);PutBit(0);PutVLIC(L);
   for (i=4:-1:1); PutBit(bitget(HL(1),i)); end;
   for l=2:L
      if (HL(l)==HL(l-1)); PutBit(0);
      elseif (HL(l)==(HL(l-1)+1)); PutBit(1);PutBit(0);
      elseif (HL(l)==(HL(l-1)-1)); PutBit(1);PutBit(1);PutBit(0);PutBit(0);
      elseif (HL(l)==(HL(l-1)+2)); PutBit(1);PutBit(1);PutBit(0);PutBit(1);
      else 
         PutBit(1);PutBit(1);PutBit(1);
         for (i=4:-1:1); PutBit(bitget(HL(l),i)); end;
      end
   end
end
if (Method==2)
   PutBit(0);PutBit(1);PutVLIC(L);
   for l=1:L
      for (i=4:-1:1); PutBit(bitget(HL(l),i)); end;
   end
end
if (Method==3)
   PutBit(1);PutBit(0);PutVLIC(L);
   for (i=4:-1:1); PutBit(bitget(HL(1),i)); end;
   r=0;
   for l=2:L
      if (HL(l)==HL(l-1)); r=r+1;
      else
         if r>0
            PutBit(0);PutBit(0);PutVLIC(r);r=0;
         end
         if (HL(l)==(HL(l-1)+1)); PutBit(0);PutBit(1); 
         elseif (HL(l)==(HL(l-1)-1)); PutBit(1);PutBit(0);
         elseif (HL(l)==(HL(l-1)+2)); PutBit(1);PutBit(1);PutBit(0);
         else 
            PutBit(1);PutBit(1);PutBit(1);
            for (i=4:-1:1); PutBit(bitget(HL(l),i)); end;
         end
      end
   end
   if r>0
      PutBit(0);PutBit(0);PutVLIC(r);r=0;
   end
end
if (Method==4)
   PutBit(1);PutBit(1);PutVLIC(L);
   for (i=5:-1:1); PutBit(bitget(HL(1),i)); end;
   r=0;
   for l=2:L
      if (HL(l)==HL(l-1)); r=r+1;
      else
         if r>0
            PutBit(0);PutVLIC(r);r=0;
         end
         if (HL(l)==(HL(l-1)+1)); PutBit(1);PutBit(0); 
         else 
            PutBit(1);PutBit(1);
            for (i=5:-1:1); PutBit(bitget(HL(l),i)); end;
         end
      end
   end
   if r>0
      PutBit(0);PutVLIC(r);r=0;
   end
end
return;  % end of PutHuffTab

function HL=GetHuffTab
global y Byte BitPos
if GetBit
   if GetBit; Method=4; else; Method=3; end;
else
   if GetBit; Method=2; else; Method=1; end;
end
L=GetVLIC;
HL=zeros(L,1);
if (Method==1)
   for (i=1:4); HL(1)=HL(1)*2+GetBit; end;
   for l=2:L
      if GetBit
         if GetBit
            if GetBit
               for (i=1:4); HL(l)=HL(l)*2+GetBit; end;
            else
               if GetBit
                  HL(l)=HL(l-1)+2;
               else
                  HL(l)=HL(l-1)-1;
               end
            end
         else
            HL(l)=HL(l-1)+1;
         end
      else
         HL(l)=HL(l-1);
      end
   end
end
if (Method==2)
   for l=1:L
      for (i=1:4); HL(l)=HL(l)*2+GetBit; end;
   end
end
if (Method==3)
   for (i=1:4); HL(1)=HL(1)*2+GetBit; end;
   l=1;
   while l<L
      l=l+1;
      if GetBit
         if GetBit
            if GetBit
               for (i=1:4); HL(l)=HL(l)*2+GetBit; end;
            else
               HL(l)=HL(l-1)+2;
            end
         else
            HL(l)=HL(l-1)-1;
         end
      else
         if GetBit
            HL(l)=HL(l-1)+1;
         else
            r=GetVLIC;
            HL(l:(l+r-1))=HL(l-1);
            l=l+r-1;
         end
      end
   end
end
if (Method==4)
   for (i=1:5); HL(1)=HL(1)*2+GetBit; end;
   l=1;
   while l<L
      l=l+1;
      if GetBit
         if GetBit
            for (i=1:5); HL(l)=HL(l)*2+GetBit; end;
         else
            HL(l)=HL(l-1)+1;
         end
      else
         r=GetVLIC;
         HL(l:(l+r-1))=HL(l-1);
         l=l+r-1;
      end
   end
end
return;  % end of GetHuffTab

% Functions to write and read a Variable Length Integer Code word
% This is a way of coding non-negative integers that uses fewer 
% bits for small integers than for large ones. The scheme is:
%   '00' +  4 bit  - integers from 0 to 15
%   '01' +  8 bit  - integers from 16 to 271
%   '10' + 12 bit  - integers from 272 to 4367
%   '11' + 16 bit  - integers from 4368 to 69903
%   not supported  - integers >= 69904
function PutVLIC(N)
global y Byte BitPos
if (N<0)
   error('huff03-PutVLIC: Number is negative.'); 
elseif (N<16)
   PutBit(0);PutBit(0);
   for (i=4:-1:1); PutBit(bitget(N,i)); end;
elseif (N<272)
   PutBit(0);PutBit(1);
   N=N-16;
   for (i=8:-1:1); PutBit(bitget(N,i)); end;
elseif (N<4368)
   PutBit(1);PutBit(0);
   N=N-272;
   for (i=12:-1:1); PutBit(bitget(N,i)); end;
elseif (N<69904)
   PutBit(1);PutBit(1);
   N=N-4368;
   for (i=16:-1:1); PutBit(bitget(N,i)); end;
else
   error('huff03-PutVLIC: Number is too large.'); 
end
return

function N=GetVLIC
global y Byte BitPos
N=0;
if GetBit
   if GetBit
      for (i=1:16); N=N*2+GetBit; end;
      N=N+4368;
   else
      for (i=1:12); N=N*2+GetBit; end;
      N=N+272;
   end
else
   if GetBit
      for (i=1:8); N=N*2+GetBit; end;
      N=N+16;
   else
      for (i=1:4); N=N*2+GetBit; end;
   end
end
return

% Functions to write and read a Bit
function PutBit(Bit)
global y Byte BitPos
BitPos=BitPos-1;
if (~BitPos); Byte=Byte+1; BitPos=8; end; 
y(Byte) = bitset(y(Byte),BitPos,Bit);
return
   
function Bit=GetBit
global y Byte BitPos
BitPos=BitPos-1;
if (~BitPos); Byte=Byte+1; BitPos=8; end; 
Bit=bitget(y(Byte),BitPos);
return;
   
