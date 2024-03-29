function [ranking, scores] = calc_scores(method, residuals)

load metadata.mat
switch method
    
    case 'none' % no weighting; a simple sum
        sum = zeros(1,numcontigs);
        for i = 1:numcontigs
            res = residuals{i};
            [~,n] = size(res);
            for j = 1:n
                sum(i) = sum(i) + norm(res(:,j));
            end
        end
        [scores, ranking] = sort(sum, 'ascend');
        
    case 'equal'
        scores = zeros(1,numcontigs);
        r1 = residuals{1};
        [~,n] = size(r1);
        for i = 1:n
            temp = zeros(1,numcontigs);
            for j = 1:numcontigs
                res = residuals{j};
                temp(j) = norm(res(:,i));
            end
            [~,idx] = sort(temp, 'ascend');
            
            for j = 1:numcontigs
                scores(idx(j)) = scores(idx(j)) + j;
            end
        end
        [scores, ranking] = sort(scores, 'ascend');
end

end