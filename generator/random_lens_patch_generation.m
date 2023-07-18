%% random lens unit generation
% parameter setup
z=zeros(1601); % total canvas size
z2=z; % copy of canvas for pixel distance calculation
z(781:820,781:820)=1; % assign first lens effective area
z2(801,801)=1; % assign first lens coordinate
d=1; % stop criteria
listx=801; % list of lens coordinates
listy=801;
idx=1; % number of lenses
%% generate lens array
while d==1
    x=round(cos(rand*pi)*randi(750)+801); % random generate lens coordinate
    y=round(cos(rand*pi)*randi(750)+801);
%     x=round(rand*900+51);
%     y=round(rand*900+51);
    if z2(x,y)==0 % calulcate nearby lens distances to check if current position is valid for new lens 
        p1=z2(x-20:x+19,y-20:y+19);
        p2=z2(x-40:x+39,y-40:y+39);
        if isempty(find(p1>0))&&(isempty(find(p2>0))==0)
            idx=idx+1
            listx=[listx,x];
            listy=[listy,y];
            z(x-19:x+19,y-19:y+19)=idx;%rand*0.9+0.1;
            z2(x,y)=1;
%             dt=z(101:900,101:900);
%             if isempty(find(dt==0))
%                 d=0;
%             end
            if (mod(idx,50)==0) % plot new lens map and update coordinates
                imagesc(z)
                pause(0.2)
            end
        end
    end
end

area=zeros(1,length(listx));
for idx=1:length(listx)
    area(idx)=length(find(z==idx));
end
            
%         p=z(x-30:x+29,y-30:y+29);
%         c1=length(find(p>0));
%         c2=length(find(p==0));
%         if (c1>0)&&(c2>600)
%             z_temp=z;
%             z_temp(x-100:x+99,y-100:y+99)=idx;
%             new_part=1;
%             for tidx=1:idx
%                 t_zone=z_temp(listx(tidx)-100:listx(tidx)+99,listy(tidx)-100:listy(tidx)+99);
%                 t_zone2=t_zone==tidx;
%                 CC=bwconncomp(t_zone2,4);
%                 max_length=0;
%                 for tidx2=1:length(CC.PixelIdxList)
%                     if length(CC.PixelIdxList{tidx2})>max_length
%                         max_length=length(CC.PixelIdxList{tidx2});
%                     end
%                 end
%                 if max_length<8000
%                     new_part=0;
%                     break;
%                 end
%             end
%             if new_part==1
%                 idx=idx+1
%                 listx=[listx,x];
%                 listy=[listy,y];
%                 z(x-100:x+99,y-100:y+99)=idx;
%                 dt=z(101:2100,101:2100);
%                 if isempty(find(dt==0))
%                     d=0;
%                 end
%                 imagesc(z)
%                 pause(0.5)
%             end
%         end
%     end
% end

% this will make lens uniformly distribute, with reasonable area for each
% lens
%you can manually add lens to spaces after the program fill in most areas
