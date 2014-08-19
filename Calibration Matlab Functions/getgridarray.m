%i starts at 1
function arrayout = getgridarray(point, gridspace, N)

arrayout = point-gridspace*(N-1)/2:gridspace:point+gridspace*(N-1)/2;

end