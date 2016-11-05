function cmap = redblue(varagin)

if nargin == 0
    n = 64;
else
    n = varagin{1};
end

vlin = linspace(0,1,n/2)'; vinv = flipud(vlin); vone = 1+0*vlin;
cmap = [vlin,vlin,vone ;
        vone,vinv,vinv ];