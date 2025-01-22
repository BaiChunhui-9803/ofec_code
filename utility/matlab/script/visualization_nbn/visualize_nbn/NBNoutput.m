function []= NBNoutput(EdgeB,EdgeBFile,EdgeA, xi,yi,zi,colorR,colorG,colorB, markerPtrs)
%EdgeBFile='EdgeB.txt';
writematrix(EdgeB,EdgeBFile);
writematrix(EdgeA,'EdgeA.txt');
writematrix(xi,'xi.txt');
writematrix(yi,'yi.txt');
writematrix(zi,'zi.txt');
writematrix(colorR,'colorR.txt');
writematrix(colorG,'colorG.txt');
writematrix(colorB,'colorB.txt');
writematrix(markerPtrs,'markerPtrs.txt');
end 