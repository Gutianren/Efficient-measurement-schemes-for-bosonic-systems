function ImportanceSampling()
    global observ len pr Nq m meas T
    
    len = m;
    pr = abs(observ(:,1)) / sum(abs(observ(:,1)));
    pr = pr';
    meas = sign(observ(:,2:Nq+1))';
    
%     for i = 1:size(meas,1)
%         for j = 1:size(meas,2)
%             if meas(i,j) == 0
%                 meas(i,j) = 1;
%             end
%         end
%     end


%     measAndP(1:Nq, 1:len) = meas(:,1:len);
%     measAndP(Nq+1, 1:len) = pr(1:len);
%     measAndP = measAndP(:,1:len);
%     measAndP = sort_meas(measAndP);%sort
%     meas = measAndP(1:Nq,1:len);
%     pr = measAndP(Nq+1,1:len);
end
