%__________________________________________________________________________
% EXTRACT onsets and RTs for first level analysis as txt, from
% presentation logs
% 
% by Roberta Bianco 2017
%__________________________________________________________________________

%===========================================================================
% 
%===========================================================================

function [stimsTime_sec] = create_onsets_RT_fromLogs
Vp = { '01', '02'};

Task = { 'visual1', 'visual2'};
Csv= {'visual2.csv', 'visual2.csv'};
dir = 'D:/Roberta/DATA/MRI/VP';

for t = 1:length(Task)
    if strcmp(Task(t), 'visual1'); bl = 'sess1'; dur = 0; Csv = 'a.csv';
    elseif strcmp(Task(t), 'visual2'); bl = 'sess2'; dur = 0; Csv = 'b.csv';
    end
    
    for Z1 = 1:length(Vp)
        disp(['----- VP_', Vp{Z1}, ' (', bl, ') -----']);
        
        % ----- Finding Nans
        RT = fopen([dir 'analysis_of_VP', Vp{Z1},'-', Csv], 'r');
        D = textscan(RT,'%*s %*s %s %s %s %s %*s %*s %*s %*s %*s %*s %*s %*s', 'delimiter', ',');
        fclose(RT);
        disp('- reading CSV');
        
        Actual_codes=[];
        stimsTime_sec=[];
        
        code=0;
        i=1;
        folder = exist([ dir , Vp{Z1},'/functional/', bl]);

        disp('- writing');
        for ln=1:size(D{1},1)
            if ~isempty(strfind(D{3}{ln},'1')) && ~isempty(strfind(D{2}{ln},'5'))
                code=1;
                Actual_codes(i) =code;
                if strcmp (bl,'sess1')
                    stimsTime_sec(i)=(str2num(D{4}{ln})); % the unit is in sec
                elseif strcmp(bl,'sess2')
                    stimsTime_sec(i)=(str2num(D{4}{ln}));
                end
                i=i+1;
            elseif ~isempty(strfind(D{3}{ln},'3')) && ~isempty(strfind(D{2}{ln},'5'))
                code=2;
                Actual_codes(i) =code;
                if strcmp (bl,'sess1')
                    stimsTime_sec(i)=(str2num(D{4}{ln})); % the unit is in sec
                elseif strcmp(bl,'sess2')
                    stimsTime_sec(i)=(str2num(D{4}{ln}));
                end
                i=i+1;
            elseif ~isempty(strfind(D{3}{ln},'2')) && ~isempty(strfind(D{2}{ln},'5'))
                code=3;
                Actual_codes(i) =code;
                if strcmp (bl,'sess1')
                    stimsTime_sec(i)=(str2num(D{4}{ln})); % the unit is in sec
                elseif strcmp(bl,'sess2')
                    stimsTime_sec(i)=(str2num(D{4}{ln}));
                end
                i=i+1;
            elseif ~isempty(strfind(D{3}{ln},'4')) && ~isempty(strfind(D{2}{ln},'5'))
                code=4;
                Actual_codes(i) =code;
                if strcmp (bl,'sess1')
                    stimsTime_sec(i)=(str2num(D{4}{ln})); % the unit is in sec
                elseif strcmp(bl,'sess2')
                    stimsTime_sec(i)=(str2num(D{4}{ln}));
                end
                i=i+1;
                
            elseif ~isempty(strfind(D{2}{ln},'2'))
                switch D{3}{ln}(1)
                    case {'1'}
                        code=5;
                    case {'3'}
                        code=6;
                end
                
                Actual_codes(i) =code;
                if strcmp (bl,'sess1')
                    stimsTime_sec(i)=(str2num(D{4}{ln})); % the unit is in sec
                elseif strcmp(bl,'sess2')
                    stimsTime_sec(i)=(str2num(D{4}{ln}));
                end
                i=i+1;
            end
            
            Aons=stimsTime_sec(Actual_codes==1);
            Cons=stimsTime_sec(Actual_codes==2);
            Bons=stimsTime_sec(Actual_codes==3);
            DOons=stimsTime_sec(Actual_codes==4);
            aons=stimsTime_sec(Actual_codes==5);
            cons=stimsTime_sec(Actual_codes==6);
            
            A= fopen([dir, Vp{Z1},'/functional/', bl, '/tcl_RT.txt'], 'wt');
            fprintf(A, '%.3d\n', Aons);

            C= fopen([dir, Vp{Z1},'/functional/', bl, '/ncl_RT.txt'], 'wt');
            fprintf(C, '%.3d\n', Cons);

            B= fopen([dir, Vp{Z1},'/functional/', bl, '/til_RT.txt'], 'wt');
            fprintf(B, '%.3d\n', Bons);

            DO= fopen([dir, Vp{Z1},'/functional/', bl, '/nil_RT.txt'], 'wt');
            fprintf(DO, '%.3d\n', DOons);
            
            a= fopen([dir, Vp{Z1},'/functional/', bl, '/tcs_RT.txt'], 'wt');
            fprintf(a, '%.3d\n', aons);

            c= fopen([dir , Vp{Z1},'/functional/', bl, '/ncs_RT.txt'], 'wt');
            fprintf(c, '%.3d\n', cons);

            fclose('all');
            
        end
        disp('--> DONE');
    end
end







