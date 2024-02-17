%%
close all;clear;
% sta_num=2;% number of states
data_filelocation = dir(['file location','\*.txt']);% txt files locations
filelocation_size = max(size(data_filelocation));% get numbers of txt files
Leakage = 0.1;% a factor for spectral cross-talk arising from donor fluorescence leakage in the acceptor channel
correctY = 1.23;% a detection correction factor of donor and acceptor channel
%%
for tracenumber = 1:1:(filelocation_size)
    data_fid = fopen([data_filelocation(tracenumber,1).folder,'\',data_filelocation(tracenumber,1).name],'r');
    F = textscan(data_fid,'%f %f %f','delimiter',',','HeaderLines',1);% read txt files
    photo_time = F{1};
    photo_donor = F{2};
    photo_acceptor = F{3};
    photo_fret = (photo_acceptor - Leakage*photo_donor)./((photo_acceptor - Leakage*photo_donor) + photo_donor/correctY);
    tempicture = figure;
    subplot(2,1,1)
    set(tempicture, 'unit', 'normalized', 'position', [0,0,1,1]);
    plot(photo_time,photo_donor,photo_time,photo_acceptor);
    xlabel('Frame');
    ylabel('A.U.');
    subplot(2,1,2);
    plot(photo_time,photo_fret);
    xlabel('Frame');
    ylabel('FRET');
    ylim([-0.2,1.4])
    zoom on;
    pause;
    zoom off;
    abc=ginput(2);
    first = ceil(abc(1,1))-photo_time(1,1);
    second = floor(abc(2,1))-photo_time(1,1);
    if second - first > 5
        writematrix(photo_fret(first:second),'intensity.txt','Delimiter','tab','WriteMode','append');
        fittingnumber = input('Whether or not HMM fitting');
        close(tempicture);
        if fittingnumber == 1
            sta_num = input('Please enter the number of states to be fitted:');
            [A,B] = HMM(photo_fret(first:second),photo_time(first:second),sta_num);
            [transitiondata] = statetransitions(A,B);
            temmatrix = [A{1,1}(:,2),B'];
            writematrix(transitiondata(2:end,:),'transitionsdata.txt','Delimiter','tab','WriteMode','append');
            writematrix(temmatrix,'HMMrawdata.txt','Delimiter','tab','WriteMode','append');
        end
    end
    %close all;
end
fclose(datatfid);