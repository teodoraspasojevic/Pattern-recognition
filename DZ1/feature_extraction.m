function [F] = feature_extraction(X)

[Xx,Xy,Xp,Xpx,Xpy] = data_extraction(X);
F(1) = max(Xpx);           % max pozicija po x osi
F(2) = min(Xpx);           % min pozicija po x osi
F(3) = max(Xpy);           % max pozicija po y osi
F(4) = min(Xpy);           % min pozicija po y osi
F(5) = F(1) - F(2);        % opseg pozicija po x osi
F(6) = F(3) - F(4);        % opseg pozicija po y osi
F(7) = max(Xp);            % max pritisak
F(8) = min(Xp);            % min pritisak
F(9) = F(7) - F(8);        % opseg pritiska
F(10) = mean(Xx);          % sr. vr. brzine po x osi
F(11) = mean(Xy);          % sr. vr. brzine po y osi

% broj odbiraka u donjoj polovini slike
F(12) = histcounts2(Xpx,Xpy,[F(2) F(1)],[F(4) (F(3)+F(4))/2]);
% broj odbiraka u gornjoj polovini slike
F(13) = histcounts2(Xpx,Xpy,[F(2) F(1)],[(F(3)+F(4))/2+1 F(3)]);
% broj odbiraka u levoj polovini slike
F(14) = histcounts2(Xpx,Xpy,[F(2) (F(1)+F(2))/2],[F(4) F(3)]);
% broj odbiraka u desnoj polovini slike
F(15) = histcounts2(Xpx,Xpy,[(F(1)+F(2))/2+1 F(1)],[F(4) F(3)]);

end