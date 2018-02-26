%clear
function rc = fileline(a,b)

filepath=a;
name=b;

%rips ./data/D1A2K D1A2K-a0a

dir=fileparts(mfilename('fullpath'));
javaaddpath([dir,'/javaplex/lib/javaplex.jar']);
import edu.stanford.math.plex4.*;
javaaddpath([dir,'/javaplex/lib/plex-viewer.jar']);
import edu.stanford.math.plex_viewer.*;
addpath([dir,'/javaplex/utility']);


formatSpec = '%d %d %f %f %f';
sizeA = [5,Inf];

%formatSpec = '%d %d %f %f %f %f';
%sizeA = [6,Inf];

PROELE = {'C','N','O','S','H'};
LIGELE = {'C','N','O','S','P','F','Cl','Br','I','H'};


fprintf('Starting Rips-Complex Calculations...... \n')

    Name=strcat(name,'.pts');
    %Name=strcat(name,'_charge.pts');
    fileID = fopen(strcat(filepath,'/pts_dir/',Name), 'r');
    B = fscanf(fileID, formatSpec, sizeA);
    fclose(fileID);
    %g=A(1,:) 
    %atom=floor(rand(1,size(A,2))*3);
    %B=[A(1,:);atom;A(2:4,:)];
    atomtypes=unique(B(2,:));
    for a0=1:length(atomtypes)
        atoms0=find(B(1,:)==0 & B(2,:)==atomtypes(a0));
        for a1=1:length(atomtypes)
            %d=[a0,a1];
            
            %[atomtypes(a0) atomtypes(a1)];
            %[B(:,atoms0(1)) B(:,atoms1(1))]
            atoms1=find(B(1,:)==1 & B(2,:)==atomtypes(a1));
            clear distances
            dist=ones(length(atoms0),length(atoms1))*100.0;
            dist0=ones(length(atoms0),length(atoms0))*100.0;
            dist1=ones(length(atoms1),length(atoms1))*100.0;

            %distances=ones(size(B,2),size(B,2))*100.0;
            %size(distances)
            for i=1:length(atoms0)
                ii=atoms0(i);
                for j=1:length(atoms1)
                    jj=atoms1(j);
                    dis = sqrt((B(3,ii) - B(3,jj))^2 + (B(4,ii) - B(4,jj))^2 + (B(5,ii) - B(5,jj))^2);
                    dist(i,j) = dis;
                end
            end
            distances=[dist0 dist;dist',dist1];

            %continue
            m_space = metric.impl.ExplicitMetricSpace(distances);
            stream = api.Plex4.createVietorisRipsStream(m_space, 1, 12, 1000);
            persistence = api.Plex4.getModularSimplicialAlgorithm(1, 2);
            intervals = persistence.computeIntervals(stream);
            %barcode_fig_outfile=['rips',int2str(a0),'_',int2str(a1),'.png'];
            %if ~exist(barcode_fig_outfile,'file')
            %        options.filename = barcode_fig_outfile;
            %        options.max_filtration_value = 12;
            %       options.max_dimension = 2;
            %       options.caption=[int2str(a0),'-',int2str(a1)];
            %       options.file_format='png';
            %       plot_barcodes(intervals, options);
            %       close
            %end
            fileID = fopen(strcat(filepath,'/bar_files/',name,'_',PROELE{a0},'-',LIGELE{a1},'.bar'), 'w');
            endpoints = homology.barcodes.BarcodeUtility.getEndpoints(intervals, 0, false);
            if size(endpoints,1) > 0

                size(endpoints,1);
                dims = zeros(1,size(endpoints,1));
                bars = [dims; endpoints(:,1)'; endpoints(:,2)'];
                fprintf(fileID, '%d %4.4f %4.4f\n', bars);
            end

            fclose(fileID);
            clear m_space;
            clear stream;
            clear persistence;
            clear intervals;
            clear endpoints;
        end
        %atomtypes(end)=[]
    end
%delete(strcat(filepath,'/bar_files/*.bar'))
    clear B
    clear atomtypes
    clear a0
    clear a1
    clear atoms0
    clear atoms1
    clear dist
    clear dist1
    clear dist0


    formatSpec = '%d %d %f %f %f %f';
    sizeA = [6,Inf];


    Name=strcat(name,'_charge.pts');
    fileID = fopen(strcat(filepath,'/pts_dir/',Name), 'r');
    B = fscanf(fileID, formatSpec, sizeA);
    fclose(fileID);

    
    atomtypes=unique(B(2,:));
    for a0=1:length(atomtypes)
        atoms0=find(B(1,:)==0 & B(2,:)==atomtypes(a0));
        for a1=1:length(atomtypes)
      
            %[atomtypes(a0) atomtypes(a1)];
            atoms1=find(B(1,:)==1 & B(2,:)==atomtypes(a1));
            clear distances
            dist=ones(length(atoms0),length(atoms1))*100.0;
            dist0=ones(length(atoms0),length(atoms0))*100.0;
            dist1=ones(length(atoms1),length(atoms1))*100.0;

            for i=1:length(atoms0)
                ii=atoms0(i);
                for j=1:length(atoms1)
                    jj=atoms1(j);
                                      
                    dis = sqrt((B(3,ii) - B(3,jj))^2 + (B(4,ii) - B(4,jj))^2 + (B(5,ii) - B(5,jj))^2);
                    
                    qi=B(6,i);
                    qj=B(6,j);
                    %dist(i,j) = dis;
                    dist(i,j) = 1.0/(1.0+exp(-100.0*qi*qj/dis));
                end
            end
            distances=[dist0 dist;dist',dist1];
            

            m_space = metric.impl.ExplicitMetricSpace(distances);
            stream = api.Plex4.createVietorisRipsStream(m_space, 1, 2.0, 1000);
            %stream = api.Plex4.createVietorisRipsStream(m_space, 1, 12, 1000);
            persistence = api.Plex4.getModularSimplicialAlgorithm(1, 2);
            intervals = persistence.computeIntervals(stream);
            %barcode_fig_outfile=['rips',int2str(a0),'_',int2str(a1),'.png'];
            %if ~exist(barcode_fig_outfile,'file')
            %        options.filename = barcode_fig_outfile;
            %        options.max_filtration_value = 12;
            %       options.max_dimension = 2;
            %       options.caption=[int2str(a0),'-',int2str(a1)];
            %       options.file_format='png';
            %       plot_barcodes(intervals, options);
            %       close
            %end
            
            fileID = fopen(strcat(filepath,'/bar_files/',name,'_',PROELE{a0},'-',LIGELE{a1},'_charge.bar'), 'w');
            endpoints = homology.barcodes.BarcodeUtility.getEndpoints(intervals, 0, false);
            if size(endpoints,1) > 0
                
                size(endpoints,1);
                dims = zeros(1,size(endpoints,1));
                bars = [dims; endpoints(:,1)'; endpoints(:,2)'];
                fprintf(fileID, '%d %4.4f %4.4f\n', bars);
            end
            fclose(fileID);
            clear m_space;
            clear stream;
            clear persistence;
            clear intervals;
            clear endpoints;
        end
        %atomtypes(end)=[]
    end
























fprintf('Done! \n')
%exit
