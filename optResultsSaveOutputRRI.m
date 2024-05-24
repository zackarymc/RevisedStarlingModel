% script optResultsSaveOutputRRI.m Saves optimized results to file
% Instigated from function revStarOptSelectedRRIfunc.m

%% File output

    allvars = whos;
    tosave = cellfun(@isempty, regexp({allvars.class}, '^matlab\.(ui|graphics)\.'));
    save(outfilename, allvars(tosave).name);%saves all variables to file

    %Append results to Excel spreadsheet

    XLfilename=strcat(xlfname,'a','.xls');%sets Excel data filename for output

    if exist(XLfilename,'file')==2 %file exists already--need to find row to append to
        XLexist=1;
        olddata=readtable(XLfilename);
        datarange=size(olddata);
        lastrow=datarange(1);%picks row dimension
        lastrowindex=strcat('A',num2str(lastrow+2));%gets the index of next free row
    else
        lastrowindex='A1';%brand new file to be written to
        XLexist=0;
    end

    if XLexist==0 %file needs to be created and data titles written
        headercells=[];
        dataidtitles1={'Analysis program' 'Dataset name' 'N datapoints' 'Ntrials'};
        dataidtitles2={'Hct(0)' 'Hct(end)'  };
        dataidtitles3={'P0final' 'VT0final' 'ChSlfinal' ...
             'XJvfinal' 'XJsfinal' 'Ppfinal' ...
             'phimaxfinal' 'PosSlfinal'...
            'net delta Vp' 'net delta Vi' 'summed Jv' 'summed Js'};
        dataidtitles4={'summed Jl' 'summed Juf' 'delta plasma protein'...
            'delta interstit protein' 'variance best'};

        headercells=[dataidtitles1 dataidtitles2 dataidtitles3 dataidtitles4];
        Theader=cell2table(headercells);
        %writetable(Theader,XLfilename,'WriteVariableNames',false,'Range','A1');
        writetable(Theader,XLfilename,'FileType','spreadsheet','WriteVariableNames',false,'Range','A1');
        lastrowindex='A2';
    end 

    %Now set up data to write to Excel
    datarow=[];
    data1={ModelID outfilename ndatatotal Ntrials};
    data2={HctEst0 HctEstlast};
    data3={P0final VT0final ChSlfinal ...
        XJvfinal XJsfinal Ppfinal ...
        phimaxfinal PosSlfinal ...
        netdeltaVpfluxes netdeltaVifluxes summedJv  summedJs ...
        summedJl summedAppliedFlux net_plasma_protein_change_fluxes ...
        net_interstitial_protein_fluxes...
        bestfittedvariance};
    datarow=[data1 data2 data3];
    Tdata=cell2table(datarow);

    writetable(Tdata,XLfilename,'FileType','spreadsheet','WriteVariableNames',false,'Range',lastrowindex);