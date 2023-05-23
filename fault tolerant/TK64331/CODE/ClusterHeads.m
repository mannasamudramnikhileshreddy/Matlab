function [CHs,rgi,rg] = ClusterHeads(gridtype,n,rnum)

r1 = 1;r2 = 1;r3 = 1;r4 = 1;r5 = 1;r6 = 1;
r7 = 1;r8 = 1;r9 = 1;r10 = 1;r11 = 1;r12 = 1;
r13 = 1;r14 = 1;r15 = 1;r16 = 1;r17 = 1;r18 = 1;

for i3 = 1:n
    if gridtype(i3) == 1
        rgi1(r1) = i3;
        rg1(r1) = rnum(i3);
        r1 = r1+1;
    elseif gridtype(i3) == 2
        rgi2(r2) = i3;
        rg2(r2) = rnum(i3);
        r2 = r2+1;
    elseif gridtype(i3) == 3
        rgi3(r3) = i3;
        rg3(r3) = rnum(i3);
        r3 = r3+1;
    elseif gridtype(i3) == 4
        rgi4(r4) = i3;
        rg4(r4) = rnum(i3);
        r4 = r4+1;
    elseif gridtype(i3) == 5
        rgi5(r5) = i3;
        rg5(r5) = rnum(i3);
        r5 = r5+1;
    elseif gridtype(i3) == 6
        rgi6(r6) = i3;
        rg6(r6) = rnum(i3);
        r6 = r6+1;
    elseif gridtype(i3) == 7
        rgi7(r7) = i3;
        rg7(r7) = rnum(i3);
        r7 = r7+1;
    elseif gridtype(i3) == 8
        rgi8(r8) = i3;
        rg8(r8) = rnum(i3);
        r8 = r8+1;
    elseif gridtype(i3) == 9
        rgi9(r9) = i3;
        rg9(r9) = rnum(i3);
        r9 = r9+1;
    elseif gridtype(i3) == 10
        rgi10(r10) = i3;
        rg10(r10) = rnum(i3);
        r10 = r10+1;
    elseif gridtype(i3) == 11
        rgi11(r11) = i3;
        rg11(r11) = rnum(i3);
        r11 = r11+1;
    elseif gridtype(i3) == 12
        rgi12(r12) = i3;
        rg12(r12) = rnum(i3);
        r12 = r12+1;
    elseif gridtype(i3) == 13
        rgi13(r13) = i3;
        rg13(r13) = rnum(i3);
        r13 = r13+1;
    elseif gridtype(i3) == 14
        rgi14(r14) = i3;
        rg14(r14) = rnum(i3);
        r14 = r14+1;
    elseif gridtype(i3) == 15
        rgi15(r15) = i3;
        rg15(r15) = rnum(i3);
        r15 = r15+1;
    elseif gridtype(i3) == 16
        rgi16(r16) = i3;
        rg16(r16) = rnum(i3);
        r16 = r16+1;
    elseif gridtype(i3) == 17
        rgi17(r17) = i3;
        rg17(r17) = rnum(i3);
        r17 = r17+1;
    elseif gridtype(i3) == 18
        rgi18(r18) = i3;
        rg18(r18) = rnum(i3);
        r18 = r18+1;
    end
end
rgi = cat(2,rgi1,rgi2,rgi3,rgi4,rgi5,rgi6...
            ,rgi7,rgi8,rgi9,rgi10,rgi11,rgi12,...
            rgi13,rgi14,rgi15,rgi16,rgi17,rgi18);
        
rg = cat(2,rg1,rg2,rg3,rg4,rg5,rg6...
            ,rg7,rg8,rg9,rg10,rg11,rg12,...
            rg13,rg14,rg15,rg16,rg17,rg18);


for i41 = 1:length(rg1)
    if min(rg1) == rg1(i41)
        rCH1 = rg1(i41);
        CH1 = rgi1(i41);
    end
end
for i42 = 1:length(rg2)
    if min(rg2) == rg2(i42)
        rCH2 = rg2(i42);
        CH2 = rgi2(i42);
    end
end
for i43 = 1:length(rg3)
    if min(rg3) == rg3(i43)
        rCH3 = rg3(i43);
        CH3 = rgi3(i43);
    end
end
for i44 = 1:length(rg4)
    if min(rg4) == rg4(i44)
        rCH4 = rg4(i44);
        CH4 = rgi4(i44);
    end
end
for i45 = 1:length(rg5)
    if min(rg5) == rg5(i45)
        rCH5 = rg5(i45);
        CH5 = rgi5(i45);
    end
end
for i46 = 1:length(rg6)
    if min(rg6) == rg6(i46)
        rCH6 = rg6(i46);
        CH6 = rgi6(i46);
    end
end
for i47 = 1:length(rg7)
    if min(rg7) == rg7(i47)
        rCH7 = rg7(i47);
        CH7 = rgi7(i47);
    end
end
for i48 = 1:length(rg8)
    if min(rg8) == rg8(i48)
        rCH8 = rg8(i48);
        CH8 = rgi8(i48);
    end
end
for i49 = 1:length(rg9)
    if min(rg9) == rg9(i49)
        rCH9 = rg9(i49);
        CH9 = rgi9(i49);
    end
end
for i410 = 1:length(rg10)
    if min(rg10) == rg10(i410)
        rCH10 = rg10(i410);
        CH10 = rgi10(i410);
    end
end
for i411 = 1:length(rg11)
    if min(rg11) == rg11(i411)
        rCH11 = rg11(i411);
        CH11 = rgi11(i411);
    end
end
for i412 = 1:length(rg12)
    if min(rg12) == rg12(i412)
        rCH12 = rg12(i412);
        CH12 = rgi2(i412);
    end
end
for i413 = 1:length(rg13)
    if min(rg13) == rg13(i413)
        rCH13 = rg13(i413);
        CH13 = rgi13(i413);
    end
end
for i414 = 1:length(rg14)
    if min(rg14) == rg14(i414)
        rCH14 = rg14(i414);
        CH14 = rgi14(i414);
    end
end
for i415 = 1:length(rg15)
    if min(rg15) == rg15(i415)
        rCH15 = rg15(i415);
        CH15 = rgi15(i415);
    end
end
for i416 = 1:length(rg16)
    if min(rg16) == rg16(i416)
        rCH16 = rg16(i416);
        CH16 = rgi16(i416);
    end
end
for i417 = 1:length(rg17)
    if min(rg17) == rg17(i417)
        rCH17 = rg17(i417);
        CH17 = rgi17(i417);
    end
end
for i418 = 1:length(rg18)
    if min(rg18) == rg18(i418)
        rCH18 = rg18(i418);
        CH18 = rgi18(i418);
    end
end

CHs = cat(2,CH1,CH2,CH3,CH4,CH5,CH6,...
            CH7,CH8,CH9,CH10,CH11,CH12,...
            CH13,CH14,CH15,CH16,CH17,CH18);

end