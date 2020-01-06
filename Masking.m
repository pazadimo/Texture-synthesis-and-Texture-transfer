function mask = Masking(E)
           E_size = size(E);           
           path = zeros(E_size(1), E_size(2));
           Energy = E;
            
           for j=2:1:E_size(1)
               for i=1:1:E_size(2)
                   if( i ==1)
                       [temp, temp_index] = min(Energy(j-1,1:2));
                       Energy(j,i) = temp + Energy(j,i);
                       path(j,i)=(temp_index - 1);
                   elseif(i ==E_size(2))
                       [temp, temp_index] = min(Energy(j-1,(E_size(2)-1):E_size(2)));
                       Energy(j,i) = temp + Energy(j,i);
                       path(j,i)=(temp_index - 2);
                   else
                       s= i-1;
                       e=i+1;
                       vectorr = Energy(j-1,s:1:e);
                       [temp, temp_index] = min(vectorr);
                       Energy(j,i) = temp + Energy(j,i);
                       path(j,i)=(temp_index - 2);
                   end
                   %temp_index
                   %i
               end
           end
           mask = zeros(E_size(1),E_size(2));
           [temp , temp_index]= min(Energy(length(Energy),:));
           for j= E_size(1):-1:2
               mask(j,(temp_index):E_size(2)) = 1;
               temp_index = path(j,temp_index)+temp_index;
           end
           mask(j,(temp_index):E_size(2)) = 1;
           
end