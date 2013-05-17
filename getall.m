l = load('all_mat_list.txt');
dead = load('21dead.txt');




parfor i = 1:425
    fid=fopen('result425.csv','a');
    isdead = 0;
    for j = 1:21
        if dead(j) == l(i)
            isdead = 1
            break;
        end
    end
	name = num2str(l(i));
    data = load(name);
    data = struct2cell(data);
    data = cell2mat(data);
    for j=1:3
    	disp(['Process ',name,' => ',num2str(j)]);
        [hrv,qtv,radiov,r,Rnum] = TEST(name,j,data);
        if ~isstruct(hrv)
        	continue;
        end
        
        tmp = sprintf('%s,%d,%d,',name,j,Rnum);
        
        tmp = strcat(tmp,print_rv(hrv),',',print_rv(qtv),',',print_rv(radiov),',');

		tmp = strcat(tmp,sprintf('%f,%f,',r.lin(1),   r.lin(2)));
		tmp = strcat(tmp,sprintf('%f,%f,',r.hyp(1),   r.hyp(2)));
		tmp = strcat(tmp,sprintf('%f,%f,',r.par(1),   r.par(2)));
		tmp = strcat(tmp,sprintf('%f,%f,',r.log(1),   r.log(2)));
		tmp = strcat(tmp,sprintf('%f,%f,',r.shlog(1), r.shlog(2)));
		tmp = strcat(tmp,sprintf('%f,%f,',r.exp(1),   r.exp(2)));
		tmp = strcat(tmp,sprintf('%f,%f,',r.atan(1),  r.atan(2)));
		tmp = strcat(tmp,sprintf('%f,%f,',r.htan(1),  r.htan(2)));
		tmp = strcat(tmp,sprintf('%f,%f,',r.ahs(1),   r.ahs(2)));
		tmp = strcat(tmp,sprintf('%f,%f,',r.ahc(1),   r.ahc(2)));
		tmp = strcat(tmp,sprintf('%d',isdead));
        
        fprintf(fid,'%s\n',tmp);
    end
    fclose(fid);
end