Xsize=length(X0)/3;

Welts =  grab_elts(node_dofs,Xsize,2,5);
SVelts =  grab_elts(node_dofs,Xsize,3,2);

Welts = repmat(Welts,1,3);
SVelts = repmat(SVelts,1,3);

WX0 = X0(Welts,Welts);
WX1 = X1(Welts,Welts);
WX2 = X2(Welts,Welts);

SVX2 = X2(SVelts,SVelts);
SVX1 = X1(SVelts,SVelts);
SVX0 = X0(SVelts,SVelts);

[Wshape,Womega] = polyeig(-WX0,-WX1,WX2);
Wevals = Womega;
Womega=Womega*alpha;
Womega = abs(Womega);
[Womega,Wix] = sort(Womega,'ascend');
Wshape = Wshape(:,Wix);
Wevals = Wevals(Wix);
Wfreqs = abs(Wevals);

[SVshape,SVomega] = polyeig(-SVX0,-SVX1,SVX2);
SVevals = SVomega;
SVomega=SVomega*alpha;
SVomega = abs(SVomega);
[SVomega,SVix] = sort(SVomega,'ascend');
SVshape = SVshape(:,SVix);
SVevals = SVevals(SVix);
SVfreqs = abs(SVevals);