function bounds = define_layers(spacing,num_channels,exp_path,save_layers)

% save_layers = 0 for N, 1 for Y 

cd(exp_path)

%L4
disp ('Look at the CSD and estimate which channels correspond with which layers')
st_4 = input('First channel in L4: ','s');
end_4 = input('Last channel in L4: ','s');
chs_4=str2double(st_4):str2double(end_4);

% L2/3
est_start_23 = str2double(st_4)-200/spacing;
end_23 = str2double(st_4)-1;
if est_start_23 <= 0
    est_start_23 = 1;
end
disp(sprintf('Estimated beginning of L2/3 so that it is <=200 microns: %d',est_start_23))
edit_23 = input('Do you want to change the start of L2/3? (Y/N): ','s');
if strcmp(edit_23,'Y')
    st_23 = input('First channel in L2/3: ','s');
    st_23 = str2double(st_23);
else
    st_23 = est_start_23;
end

%L5
st_5a = str2double(end_4)+1;
est_end_5a = str2double(end_4)+125/spacing;
if est_end_5a > num_channels            
    est_end_5a = num_channels;
end
disp(sprintf('Estimated end of L5a so that it is <=125 microns: %d',est_end_5a))
edit_5a = input('Do you want to change the end of L5a? (Y/N): ','s');
if strcmp(edit_5a,'Y')
    end_5a = input('Last channel in L5a: ','s');
    end_5a = str2double(end_5a);
else
    end_5a = est_end_5a;
end

st_5b = end_5a+1;
est_end_5b = end_5a+125/spacing;
if est_end_5b > num_channels            
    est_end_5b = num_channels;
end
disp(sprintf('Estimated end of L5b so that it is <=125 microns: %d',est_end_5b))
edit_5b = input('Do you want to change the end of L5b? (Y/N): ','s');
if strcmp(edit_5b,'Y')
    end_5b = input('Last channel in L5: ','s');
    end_5b = str2double(end_5b);
else
    end_5b = est_end_5b;
end

%L6
st_6 = end_5b+1;
if st_6 <= num_channels
    end_6 = num_channels;     % for now, ignores white matter
else
    st_6 = [];
    end_6 = [];
end

bounds = [st_23 mean([end_23 chs_4(1)]) mean([chs_4(end) st_5a]) mean([end_5a st_5b]) mean([end_5b st_6])];
bounds = bounds(~isnan(bounds));

layers = zeros(1,num_channels);
layers(st_23:end_23) = 2.5;
layers(chs_4) = 4;
layers(st_5a:end_5a) = 5;
layers(st_5b:end_5b) = 5.5;
layers(st_6:end_6) = 6;

if save_layers
    % save the layer information is relevant directories
    cd ..
    expdir = cd;
    dirname= ...
        uigetdir(expdir,'In which experiment folders do you want to save the layer information?');
        save(sprintf('%s\\layers.mat',dirname),'layers')
    moredirs = input('Do you want to save the layers in more places? Y or N: ','s');
    if strcmp(moredirs,'Y')
        howmany = input('How many more directories do you want to save the layer info in? ','s');
        for n = 1:str2double(howmany)
            dirname = uigetdir(expdir,'In which experiment folders do you want to save the layer information?');
            save(sprintf('%s\\layers.mat',dirname),'layers')
        end
    end
    cd(exp_path)
end

return
