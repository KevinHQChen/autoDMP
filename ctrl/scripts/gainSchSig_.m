function sig = gainSchSig_(ywt)
persistent sigPrev;
if isempty(sigPrev)
    sigPrev = 1;
end

if ( ywt(1) > 0.5 && all(ywt(2:end) < 0.5) )...
        || ( ywt(4) > 0.5 && all(ywt([1:3,5:end]) < 0.5) )
    sig = 1;
elseif ( ywt(2) > 0.5 && ywt(3) > 0.5 && all(ywt([1,4:end]) < 0.5) )...
        || ( ywt(5) > 0.5 && ywt(6) > 0.5 && all(ywt([1:4]) < 0.5) )
    sig = 2;
elseif ( ywt(1) > 0.5 && ywt(2) > 0.5 && all(ywt([3:end]) < 0.5) )...
        || ( ywt(4) > 0.5 && ywt(5) > 0.5 && all(ywt([1:3,6]) < 0.5) )
    sig = 3;
else
    sig = sigPrev;
end
sigPrev = sig;
