function serial = zigzag(input)

ind = reshape(1:numel(input), size(input));   
ind = fliplr( spdiags( fliplr(ind) ) );    
ind(:,1:2:end) = flipud( ind(:,1:2:end) );  
ind(ind==0) = [];   

serial = input(ind)  ;

end
