% function stop = store_obj_history(~, optimValues, state)
%     persistent fvals
%     stop = false;
% 
%     switch state
%         case 'init'
%             fvals = [];
%         case 'iter'
%             fvals(end+1) = optimValues.fval;
%         case 'done'
%             assignin('base', 'last_fval_history', fvals);
%     end
% end
