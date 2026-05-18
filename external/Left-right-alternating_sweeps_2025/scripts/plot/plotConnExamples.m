function plotConnExamples(recs)
%PLOTEXAMPLES Summary of this function goes here
%%   Detailed explanation goes here
S = SweepsSettings;

recNames.idid = ["28304_1", "28304_1", "28258_4"];
recNames.idconj = ["26035_1", "28304_1", "26648_1"];
recNames.conjgrid= ["26035_1", "28304_1", "25843_1"];
% 
uids.idid =["2_1823", "2_1403"; "2_1947", "2_0295"; "2_0937", "2_0175"];
uids.idconj =["2_0468", "2_0206"; "2_1807", "2_0807"; "1_1213", "1_1053"];
uids.conjgrid =["2_0570","2_0147"; "2_1382", "2_0584";"2_1547","2_1323"];
% 
tl = tiledlayout(3, 11, 'TileSpacing','tight', 'TileIndexing','rowmajor');
% for i = 1:3
for comb = fieldnamesstr(recNames)
    for i = 1:3
        rec = recs([recs.recName]==recNames.(comb)(i));
        units = rec.units;
        ids = uids.(comb)(i, :);
        ipair = find(rec.unitIds(:, 1)==ids(1) & rec.unitIds(:, 2)==ids(2));
        xcorrs = rec.xcorrs.raw;
        xcorrs = xcorrs(:, rec.isconn);
        xcorr = xcorrs(:, ipair);
        us(1) = units([units.id]==ids(1));
        us(2) = units([units.id]==ids(2));
        plotOnePair(us, xcorr)
        if i<3
            ax = nexttile;
            ax.Visible="off";
        end
        if i==1
        end
    end
end
end

