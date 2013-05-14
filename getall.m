l = load('all_mat_list.txt');
dead = load('21dead.txt');

fid=fopen('result426.txt','w');

set(gcf,'visible','off');

parfor i = 1:425
    isdead = 0;
    for j = 1:21
        if dead(j) == l(i)
            isdead = 1
        end
    end
	name = num2str(l(i));
    for j=1:3
    	disp(['Process ',name,' => ',num2str(j)]);
        [hrv,qtv,r] = TEST(name,j);
        if ~isstruct(r)
        	continue;
        end
        
        fprintf(fid,'%s,%i,',name,j);

        tmp = print_rv(hrv);
        fprintf(fid,'%s,',tmp);

        tmp = print_rv(qtv);
        fprintf(fid,'%s,',tmp);

		fprintf(fid,'%f,%f,',r.lin(1),   r.lin(2));
		fprintf(fid,'%f,%f,',r.hyp(1),   r.hyp(2));
		fprintf(fid,'%f,%f,',r.par(1),   r.par(2));
		fprintf(fid,'%f,%f,',r.log(1),   r.log(2));
		fprintf(fid,'%f,%f,',r.shlog(1), r.shlog(2));
		fprintf(fid,'%f,%f,',r.exp(1),   r.exp(2));
		fprintf(fid,'%f,%f,',r.atan(1),  r.atan(2));
		fprintf(fid,'%f,%f,',r.htan(1),  r.htan(2));
		fprintf(fid,'%f,%f,',r.ahs(1),   r.ahs(2));
		fprintf(fid,'%f,%f,',r.ahc(1),   r.ahc(2));
        fprintf(fid,'%d\n',isdead);
    end
end

set(figure(1),'visible','on');
set(figure(2),'visible','on');