%% code for generating FreeCAD FCMACRO for placing CAD units into array
acoordinate;
for i=1:1
%     fileID = fopen(['design1_600umw_15mmr',num2str(i),'.txt'],'w');
    fileID = fopen('design1.txt','w');
    fprintf(fileID, 'Exponential Function\n\n');
%     for idx=100*(i-1)+1:100*i
    for idx=1:2:length(a)
        fprintf(fileID, 'obj=doc.copyObject(sphere)\n');
    %     fprintf(fileID, ['obj.Radius=',num2str(rand/12.5+1.96),'\n']);
    %     fprintf(fileID, 'obj.Radius=2\n');
        fprintf(fileID, ['obj.Placement.Base=FreeCAD.Vector(',num2str(a(idx,1)),',',num2str(a(idx,2)),',0)\n']);
    end
    fclose(fileID);
end