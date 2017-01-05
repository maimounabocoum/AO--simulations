for k=1:1:99
    close all
    for j=10*k:1:10*k+10%size(phantom_wave,3)
        v(:) = signal(:,j);
        [N,M]=size(v);
        v=v/max(max(v));
        t = Temps(j,1);
        figure(j)
        %figure
        [N,M]=size(v);
        v=v/max(max(v));
        for i=1:N_elements
            plot(2*(0:N-1)/fs+t,v(:,i)+i), hold on
        end
        hold off
        title('Individual traces')
        xlabel('Time [s]')
        ylabel('Normalized response')
        axis([t t+2*N/fs 0 M+1])
        %     v(:) = signal(:,j);
        %     [N,M]=size(v);
        %     v=v/max(max(v));
        %     t = Temps(j,1);
        %     figure(j)
        %     if size(v,2)>1
        %         plot(2*(0:N-1)/fs+t,sum(v'))
        %     else
        %         plot(2*(0:N-1)/fs+t,v(:,i)+i)
        %     end
        %     title('Summed response')
        %     xlabel('Time [s]')
        %     ylabel('Normalized response')
        %     axis([t t+2*N/fs 0 M+1])
        
        w = v';
        %save(['Registered_Data'],'w');
        %ScatterInfo.fs = fs;
        %ScatterInfo.t  = t;
        disp('Data saved')
    end
        pause;
end