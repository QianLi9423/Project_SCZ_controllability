function cohi_mat_fast(readFileName,saveFileName,TR,bandPass)
X = importdata(readFileName);
cohiMatrix = wtcMatrix(X,TR,bandPass);
save(saveFileName,'cohiMatrix');
end