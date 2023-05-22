function outp = test(flex1inp,flex2inp,EMGinp)
if flex1inp >= 890 && flex1inp <= 895
    if flex2inp >= 929 && flex2inp <= 935
        if EMGinp == 16
            outp = 'Relax position';
        end
    end
elseif flex1inp >= 980 && flex1inp <= 995
    if flex2inp >= 930 && flex2inp <= 945
        if EMGinp == 20
            outp = 'Thumbs_up';
        end
    end
elseif flex1inp >= 970 && flex1inp <= 985
    if flex2inp >= 965 && flex2inp <= 985
        if EMGinp == 32
            outp = 'Fist position';
        end
    end
end
end