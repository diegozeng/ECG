% @wave => wave to be filtered
% @Rnum => RR number
% @w    => window size
function[wave] = MeanFilter(wave,Rnum,w,level)
    %w = 20;
    %level = 0.2
    avg = zeros(Rnum - 2 * w,1);
    for i = 1:w
        avg = avg + wave(i:Rnum + i - 2 * w - 1);
    end
    for i = w + 2:2 * w + 1
        avg = avg + wave(i:Rnum + i - 2 * w - 1);
    end
    w = 2 * w;
    avg = avg / w;
    up = 1 + level;
    down = 1 - level;
    for i = w:Rnum - w
        a = avg(i - w + 1);
        if(wave(i) > a * up || wave(i) < a * down)
            wave(i) = a;
        end
    end

end