function val = qsdpmlobj(Y,NY,C,b2)
%Function takes a matrix X and scalar b and outputs the value of the
%following function:
%\fnf(\vx)=\sum_{i}( \sum_{j}\vx_{i,j} \log( \sum_{j}\vx_{i,j}) ) - 
% (\sum_{i}\sum_{j}\vC_{i,j}\vx_{i,j} +\sum_{i}\sum_{j}\vx_{i,j} \log \vx_{i,j})
% support = find(Y~=0);
% Rowsum=sum(Y,2);
% Rowsumsupport= find(Rowsum~=0);
% val=sum(sum(Y(support).*(-C(support)-log(Y(support)))))+sum(Rowsum(Rowsumsupport).*log(Rowsum(Rowsumsupport)))

%val=sum(sum(Y+logerror,2).*log(sum(Y+logerror,2)))-sum(sum((Y+logerror).*log(Y+logerror)))-sum(sum((Y+logerror).*C))

%val=sum(sum(Y,2).*log(sum(Y,2)))-sum(sum((Y).*log(Y)))-sum(sum((Y).*C))
val=sum(sum(rel_entr_quad([Y,NY],[Y,NY]*ones(b2+1,b2+1))))+sum(sum(C.*Y));
%val=sum(sum(rel_entr(Y,sum(Y,2)*ones(1,b2))))-sum(sum(C.*Y))
