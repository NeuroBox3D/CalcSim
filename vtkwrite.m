function vtkwrite(x,y,z,W,time,outfolder)
        filename = (sprintf('%s/step_t%d.vtk',outfolder,time));        
        fid = fopen(filename, 'w'); 
        nr_of_elements = length(x);
        
        %ASCII file header
        fprintf(fid, '# vtk DataFile Version 3.0\n');
        fprintf(fid, 'VTK from Matlab\n');
        fprintf(fid, 'BINARY\n\n');
        fprintf(fid, 'DATASET STRUCTURED_GRID\n');
        fprintf(fid, ['DIMENSIONS ' num2str(size(x,1)) ' ' num2str(size(y,2)) ' ' num2str(size(z,3)) '\n']);
        fprintf(fid, ['POINTS ' num2str(nr_of_elements) ' float\n']);
        fclose(fid);

        %append binary x,y,z data
        fid = fopen(filename, 'a'); 
        fwrite(fid, [reshape(x,1,nr_of_elements);  reshape(y,1,nr_of_elements); reshape(z,1,nr_of_elements)],'float','b');

        %append another ASCII sub header
        fprintf(fid, ['\nPOINT_DATA ' num2str(nr_of_elements) '\n']);
        %append some scalar data
        fprintf(fid, '\nSCALARS CYT float\n'); %ASCII header
        fprintf(fid, 'LOOKUP_TABLE default\n'); %ASCII header
        fwrite (fid, reshape(W(:,1),1,nr_of_elements),'float','b'); %binary data
        
        %append some scalar data
        fprintf(fid, '\nSCALARS ER float\n'); %ASCII header
        fprintf(fid, 'LOOKUP_TABLE default\n'); %ASCII header
        fwrite (fid, reshape(W(:,2),1,nr_of_elements),'float','b'); %binary data
                
        %append some scalar data
        fprintf(fid, '\nSCALARS CALBIN float\n'); %ASCII header
        fprintf(fid, 'LOOKUP_TABLE default\n'); %ASCII header
        fwrite (fid, reshape(W(:,3),1,nr_of_elements),'float','b'); %binary data
        
        %append some scalar data
        fprintf(fid, '\nSCALARS ERBUFF float\n'); %ASCII header
        fprintf(fid, 'LOOKUP_TABLE default\n'); %ASCII header
        fwrite (fid, reshape(W(:,4),1,nr_of_elements),'float','b'); %binary data
        
        %append some scalar data
        fprintf(fid, '\nSCALARS VOLTAGE float\n'); %ASCII header
        fprintf(fid, 'LOOKUP_TABLE default\n'); %ASCII header
        fwrite (fid, reshape(W(:,5),1,nr_of_elements),'float','b'); %binary data
        
        fclose(fid);
end