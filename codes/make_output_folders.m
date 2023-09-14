function make_output_folders()

lst_folders = {'spacedata'};

for i=1:length(lst_folders)
    if ~isfolder(lst_folders{i})
        fprintf(sprintf('Making %s output folder\n',lst_folders{i}))
        mkdir(lst_folders{i})
    end
end

end