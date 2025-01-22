a=[1 2 3 4 5 6 7 8 9];
fid=fopen('test.txt','wb');%写入'w'
fwitre(fid,a,'double');
fclose(fid);

fidin=fopen('test.txt','rb')%读取'r'
[A,COUNT]=fread(fidin,'double')
fclose(fidin)