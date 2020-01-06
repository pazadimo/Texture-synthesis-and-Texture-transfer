
image2=imread('weave.jpg');
image3=im2double(image2);
Size=size(image3);
Batch_size = 100;
step = 1;
Overlap_size = floor(Batch_size/6);
h_image = Size(1);
w_image = Size(2); 
if(ndims(image3)==2)
    image=ones(Size(1),Size(2),1);
    output_image= double(ones(h_image*5, w_image*5,1));
    potential_blocks=ones((h_image - Batch_size)*(w_image - Batch_size),Batch_size,Batch_size,1);
    %potential_blocks=ones((h_image - Batch_size)*(w_image - Batch_size),Batch_size,Batch_size,1);

    image(:,:,1) = image3(:,:)
else
    output_image= double(ones(h_image*5, w_image*5,3));
    potential_blocks=ones((h_image - Batch_size)*(w_image - Batch_size),Batch_size,Batch_size,3);
    %potential_blocks=ones((h_image - Batch_size)*(w_image - Batch_size),Batch_size,Batch_size,3);

    image = image3;
end
h_N_block =  floor((h_image*5 - Overlap_size) / (Batch_size - Overlap_size));
w_N_block =  floor((w_image*5 - Overlap_size) / (Batch_size - Overlap_size));
tolerance = 0.3;
%choosing the first block randomly
first_block_h=randi(h_image - Batch_size);
first_block_w=randi(w_image - Batch_size);
output_image(1:Batch_size, 1:Batch_size, :) = image(first_block_h:(first_block_h+Batch_size-1),first_block_w:(first_block_w+Batch_size-1),:);
unoverlap = Batch_size - Overlap_size;

%potential_blocks=ones((h_image - Batch_size)*(w_image - Batch_size),Batch_size,Batch_size,3);
c=1;
indexes =0;
for j = 1:1:h_image - Batch_size % +1
               for i = 1:1:w_image - Batch_size % +1
                   index = (w_image - Batch_size)*(j-1) + i;
                   potential_blocks(index, : , :,:) = image(j:(j+Batch_size-1),i:(i+Batch_size-1),:);
               end
end
n_potential = (h_image - Batch_size)*(w_image - Batch_size);
errors = zeros(1,(h_image - Batch_size)*(w_image - Batch_size));

for j_block = 1:1:h_N_block
    j_block
    for i_block = 1:1:w_N_block
        if( j_block == 1 && i_block == 1)
            j_current = 1;
            i_current = 1;
            
        elseif( j_block == 1 && i_block < w_N_block)
            j_current = 1;
            i_current = ((i_block-1)* (Batch_size)+1);
            %choosed_block_horizon = output_image(j_current:(j_current + Batch_size -1), i_current:(i_current + Overlap_size-1),:);
        elseif( i_block == 1 && j_block < h_N_block)
            j_current = ((j_block-1)* (Batch_size)+1);
            i_current = 1;
            %choosed_block_vertical = output_image(j_current:(j_current + Overlap_size -1), i_current:(i_current + Batch_size-1),:);
        else
            j_current = ((j_block-1)* (Batch_size) +1);
            i_current = ((i_block-1)* (Batch_size) +1);

        end
        if (i_block == 1 && j_block == 1)
            output_image(1:Batch_size, 1:Batch_size, :) = output_image(1:Batch_size, 1:Batch_size, :);
        elseif ( j_block == 1 ) 
           block_index_temp = randi(length(potential_blocks));
           block_index = (block_index_temp);
           output_image(j_current:(j_current+ Batch_size-1),i_current:(i_current+ Batch_size-1),:)= potential_blocks(block_index,:,:,:);
        elseif ( i_block == 1 )
           block_index_temp = randi(length(potential_blocks));
           block_index = (block_index_temp);
           output_image(j_current:(j_current+ Batch_size-1),i_current:(i_current+ Batch_size-1),:)= potential_blocks(block_index,:,:,:);
        else 
            block_index_temp = randi(length(potential_blocks));
           block_index = (block_index_temp);
           output_image(j_current:(j_current+ Batch_size-1),i_current:(i_current+ Batch_size-1),:)= potential_blocks(block_index,:,:,:);
        end
    end
end


figure, imshow(output_image)
title("Random3 --- Patch Size = "+Batch_size)