function [state2] = statetransitions(FRET_trace,seq_show)
b = 0;
state = zeros(1,3);
for i = 2:length(seq_show)
    if seq_show(1,i) == seq_show(1,i-1)
    else
        b = i-1-b;
        state = [state;mean(FRET_trace{1,1}(i-b:i-1,2)),mean(seq_show(1,i-b:i-1)),b];
        b = i-1;
    end
end
state2 = zeros(1,5);
for i = 2:1:length(state(:,1))-1
    state2 = [state2;state(i,1),state(i,2),state(i+1,1),state(i+1,2),state(i,3)];
end

