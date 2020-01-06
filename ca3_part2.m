
image2=imread('rice.bmp');
image=im2double(image2);
image3=imread('paper.png');
image_back=im2double(image3);
Size=size(image);
Size_back=size(image_back);
Batch_size = 15;
Overlap_size =4 %floor(Batch_size/6);
h_image = Size(1);
w_image = Size(2); 
h_image_back = Size_back(1);
w_image_back = Size_back(2); 
output_image= double(ones(h_image_back, w_image_back,3));
h_N_block =  floor((h_image_back- Overlap_size) / (Batch_size - Overlap_size));
w_N_block =  floor((w_image_back - Overlap_size) / (Batch_size - Overlap_size));
tolerance = 0.3;
alpha = 0.6;
%choosing the first block randomly
first_block_h=randi(h_image - Batch_size);
first_block_w=randi(w_image - Batch_size);
output_image(1:Batch_size, 1:Batch_size, :) = image(first_block_h:(first_block_h+Batch_size-1),first_block_w:(first_block_w+Batch_size-1),:);
unoverlap = Batch_size - Overlap_size;

potential_blocks=ones((h_image - Batch_size)*(w_image - Batch_size),Batch_size,Batch_size,3);
c=1;
indexes =0;
intensity_texture = rgb2gray(image);
intensity_back = rgb2gray(image_back);

for j = 1:1:h_image - Batch_size % +1
   for i = 1:1:w_image - Batch_size % +1
       index = (w_image - Batch_size)*(j-1) + i;
       potential_blocks(index, : , :,:) = image(j:(j+Batch_size-1),i:(i+Batch_size-1),:);
   end
end
n_potential = (h_image - Batch_size)*(w_image - Batch_size);
errors = ones(1,(h_image - Batch_size)*(w_image - Batch_size))*10000;
j_i_errors = zeros(length(errors), 2);

for j_block = 1:1:h_N_block
    j_block
    for i_block = 1:1:w_N_block
        if( j_block == 1 && i_block == 1)
            j_current = 1;
            i_current = 1;
            
        elseif( j_block == 1 && i_block < w_N_block)
            j_current = 1;
            i_current = ((i_block-1)* (Batch_size - Overlap_size)+1);
            choosed_block_horizon = output_image(j_current:(j_current + Batch_size -1), i_current:(i_current + Overlap_size-1),:);
        elseif( i_block == 1 && j_block < h_N_block)
            j_current = ((j_block-1)* (Batch_size - Overlap_size)+1);
            i_current = 1;
            choosed_block_vertical = output_image(j_current:(j_current + Overlap_size -1), i_current:(i_current + Batch_size-1),:);
        else
            j_current = ((j_block-1)* (Batch_size - Overlap_size) +1);
            i_current = ((i_block-1)* (Batch_size - Overlap_size) +1);
            if((i_current + Overlap_size-1) < (1+ 5*w_image) && (j_current + Overlap_size-1) < (1+ 5*h_image) )
                choosed_block_horizon = output_image(j_current:(j_current + Batch_size -1), i_current:(i_current + Overlap_size-1),:);
                choosed_block_vertical = output_image(j_current:(j_current + Overlap_size -1), i_current:(i_current + Batch_size-1),:);
                choosed_block_over = output_image(j_current:(j_current + Overlap_size -1), i_current:(i_current + Overlap_size-1),:);
        
            end
        end
        intensity_patch_back =intensity_back(j_current:(j_current + Batch_size -1),i_current:(i_current + Batch_size -1));
        if (i_block == 1 && j_block == 1)
            output_image(1:Batch_size, 1:Batch_size, :) = output_image(1:Batch_size, 1:Batch_size, :);
        elseif ( j_block == 1 ) 
           for j = 1:10:h_image - Batch_size % +1
               for i = 1:10:w_image - Batch_size % +1
                   patch = image(j:(j+Batch_size-1),i:(i+Batch_size-1),:);
                   intensity_patch = intensity_texture(j:(j+Batch_size-1),i:(i+Batch_size-1),:);
                   index = (w_image - Batch_size)*(j-1) + i;
                   errors(1, index) = alpha *sum(sum(sum((choosed_block_horizon - (patch(1:Batch_size ,1:Overlap_size ,:))).^2)));
                   errors(1, index) = errors(1, index)+(1-alpha) *(sum(sum((intensity_patch_back(:,:) - (intensity_patch(:,:))).^2)));
                   
                   j_i_errors(index,1) = j;
                   j_i_errors(index,2) = i;
               end
           end  
           %[min_errors,p] = min(errors);
           %errors(1,p) = 10000;
           [min_errors,p] = min(errors);
           indexes = find(errors <= min_errors * (1+tolerance));
           block_index_temp = randi(length(indexes));
           block_index = indexes(block_index_temp);
           
           %%%%%%
           j_error= j_i_errors(block_index,1);
           i_error= j_i_errors(block_index,2);
           choosed_patch = image(j_error:(j_error+Batch_size-1) ,i_error:(i_error+Batch_size-1) ,:);
           E = sum((choosed_block_horizon - image(j_error:(j_error+Batch_size-1) ,i_error:(i_error+Overlap_size-1) ,:)).^2,3);
%            E_size = size(E);           
%            path = zeros(E_size(1), E_size(2));
%            Energy = E;
%             
%            for j=2:1:E_size(1)
%                for i=1:1:E_size(2)
%                    if( i ==1)
%                        [temp, temp_index] = min(Energy(j-1,1:2));
%                        Energy(j,i) = temp + Energy(j,i);
%                        path(j,i)=(temp_index - 1);
%                    elseif(i ==E_size(2))
%                        [temp, temp_index] = min(Energy(j-1,(E_size(2)-1):E_size(2)));
%                        Energy(j,i) = temp + Energy(j,i);
%                        path(j,i)=(temp_index - 2);
%                    else
%                        s= i-1;
%                        e=i+1;
%                        vectorr = Energy(j-1,s:1:e);
%                        [temp, temp_index] = min(vectorr);
%                        Energy(j,i) = temp + Energy(j,i);
%                        path(j,i)=(temp_index - 2);
%                    end
%                    temp_index
%                    i
%                end
%            end
%            mask = zeros(E_size(1),E_size(2));
%            [temp , temp_index]= min(Energy(length(Energy),:));
%            for j= E_size(1):-1:2
%                mask(j,(temp_index):E_size(2)) = 1;
%                temp_index = path(j,temp_index)+temp_index;
%            end
%            mask(j,(temp_index):E_size(2)) = 1;
           mask = Masking(E);
           Filter=ones(Batch_size,Batch_size,3);
           Filter(:,1:Overlap_size,1)=mask;
           Filter(:,1:Overlap_size,2)=mask;
           Filter(:,1:Overlap_size,3)=mask;
           %%%%%%
           
           
           %output
           x=output_image(j_current:(j_current+ Batch_size-1),i_current:(i_current+ Batch_size-1),:);
           output_image(j_current:(j_current+ Batch_size-1),i_current:(i_current+ Batch_size-1),:)=x.*(Filter==0)+ choosed_patch.*(Filter==1);
        
        elseif ( i_block == 1 )
           for j = 1:10:h_image - Batch_size % +1
               for i = 1:10:w_image - Batch_size % +1
                   index = (w_image - Batch_size)*(j-1) + i;
                   patch = image(j:(j+Batch_size-1),i:(i+Batch_size-1),:);
                   intensity_patch = intensity_texture(j:(j+Batch_size-1),i:(i+Batch_size-1),:);
                   errors(1, index) = sum(sum(sum((choosed_block_vertical - (patch(1:Overlap_size, : ,:))).^2)));
                   errors(1, index) = errors(1, index)+(1-alpha) *(sum(sum((intensity_patch_back(:,:) - (intensity_patch(:,:))).^2)));
                  
                   j_i_errors(index,1) = j;
                   j_i_errors(index,2) = i;
               end
           end 
%            [min_errors,p] = min(errors);
%            errors(1,p) = 10000;
           [min_errors,p] = min(errors);
           indexes = find(errors <= min_errors * (1+tolerance));
           block_index_temp = randi(length(indexes));
           block_index = indexes(block_index_temp);
           
           
           j_error= j_i_errors(block_index,1);
           i_error= j_i_errors(block_index,2);
           E = sum((choosed_block_vertical - image(j_error:(j_error+Overlap_size-1) ,i_error:(i_error+Batch_size-1) ,:)).^2,3);
           choosed_patch = image(j_error:(j_error+Batch_size-1) ,i_error:(i_error+Batch_size-1) ,:);

           mask_temp = Masking(E');
           mask = mask_temp';
           Filter=ones(Batch_size,Batch_size,3);
           Filter(1:Overlap_size,:,1)=mask;
           Filter(1:Overlap_size,:,2)=mask;
           Filter(1:Overlap_size,:,3)=mask;
           %output_image(j_current:(j_current+ Batch_size-1),i_current:(i_current+ Batch_size-1),:)= potential_blocks(block_index,:,:,:);
           x=output_image(j_current:(j_current+ Batch_size-1),i_current:(i_current+ Batch_size-1),:);
           output_image(j_current:(j_current+ Batch_size-1),i_current:(i_current+ Batch_size-1),:)=x.*(Filter==0)+ choosed_patch.*(Filter==1);

        else 
            
           for j = 1:10:h_image - Batch_size % +1
               for i = 1:10:w_image - Batch_size % +1
                   patch = image(j:(j+Batch_size-1),i:(i+Batch_size-1),:);
                   intensity_patch = intensity_texture(j:(j+Batch_size-1),i:(i+Batch_size-1),:);
                   
                   index = (w_image - Batch_size)*(j-1) + i;
                   errors(1, index) = alpha*sum(sum(sum((choosed_block_horizon - (patch(1:Batch_size ,1:Overlap_size ,:))).^2)));
                   errors(1, index) = errors(1, index)+(1-alpha) *(sum(sum((intensity_patch_back(:,:) - (intensity_patch(:,:))).^2)));
                   
                   errors(1, index) = errors(1, index)+ alpha*sum(sum(sum((choosed_block_vertical - (patch(1:Overlap_size, : ,:))).^2)));
                   
                   errors(1, index) = errors(1, index)- alpha*sum(sum(sum((choosed_block_over - (patch(1:Overlap_size, 1:Overlap_size ,:))).^2)));
                   j_i_errors(index,1) = j;
                   j_i_errors(index,2) = i;
               end
           end
%            [min_errors,p] = min(errors);
%            errors(1,p) = 10000;
           [min_errors,p] = min(errors);
           indexes = find(errors <= min_errors * (1+tolerance));
           block_index_temp = randi(length(indexes));
           block_index = indexes(block_index_temp);
           
           
           j_error= j_i_errors(block_index,1);
           i_error= j_i_errors(block_index,2);
           E_horizon = sum((choosed_block_horizon - image(j_error:(j_error+Batch_size-1) ,i_error:(i_error+Overlap_size-1) ,:)).^2,3);
           E_vertical = sum((choosed_block_vertical - image(j_error:(j_error+Overlap_size-1) ,i_error:(i_error+Batch_size-1) ,:)).^2,3);
           %choosed_patch = image(j_error:(j_error+Batch_size-1) ,i_error:(i_error+Batch_size-1) ,:);
           
           mask = Masking(E_horizon);
           Filter=ones(Batch_size,Batch_size,3);
           Filter(:,1:Overlap_size,1)=mask;
           Filter(:,1:Overlap_size,2)=mask;
           Filter(:,1:Overlap_size,3)=mask;
           mask_temp = Masking(E_vertical');
           mask2=mask_temp';
           Filter(1:Overlap_size,:,1)=Filter(1:Overlap_size,:,1).*mask2;
           Filter(1:Overlap_size,:,2)=Filter(1:Overlap_size,:,2).*mask2;
           Filter(1:Overlap_size,:,3)=Filter(1:Overlap_size,:,3).*mask2;
           %output_image(j_current:(j_current+ Batch_size-1),i_current:(i_current+ Batch_size-1),:)= potential_blocks(block_index,:,:,:);
           %x=output_image(j_current:(j_current+ Batch_size-1),i_current:(i_current+ Batch_size-1),:);
           output_image(j_current:(j_current+ Batch_size-1),i_current:(i_current+ Batch_size-1),:)=output_image(j_current:(j_current+ Batch_size-1),i_current:(i_current+ Batch_size-1),:).*(Filter==0)+ image(j_error:(j_error+Batch_size-1) ,i_error:(i_error+Batch_size-1) ,:).*(Filter==1);

        end
    end
end