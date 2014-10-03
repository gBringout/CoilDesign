function [result] = CalculateInductance(filename)


ax=actxserver('FastHenry2.Document');
ax.invoke('Run',['"' filename '"']);
while(ax.invoke('IsRunning'))
    pause(0.1);
end
result = ax.invoke('GetInductance');
ax.invoke('Quit');
ax=[];