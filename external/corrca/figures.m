% This code generates all figures in the paper:
%
% Lucas C. Parra, Stefan Haufe, Jacek P. Dmochowski, "Correlated Components
% Analysis --- Extracting Reliable Dimensions in Multivariate Data",
% arXiv:1801.08881, Jan 26, 2018. https://arxiv.org/abs/1801.08881
%
% Some of them take a while to generate, in particular the 'simulations'.
% So one can select just a subset or individual figures by choosing the
% desired number(s) in the loop over figure_nr
%
% Warning! The code in this file is not pretty! We support open sharing of
% code, but alas, do not have time to do it nicely. However, the core
% functions that are called here (corrca.m mcca.m lda.m, etc.) have been
% written with some care.
%
% January 19, 2019, Lucas, Jacek, Stefan

clear all
figurepath = 'figures/'; % this is where we will save figures (except, simulations, which are saved in $figurepath/simulations)
addpath('my-matlab-tools'); % for sublable.m fdr.m topoplot.m ploterp.m
addpath('stefan'); % for est_K_*.m functions

all_tests{1}='CorrCA vs FLD vs JD';
all_tests{2}='CorrCA version speed';
all_tests{3}='find shared dimension in 3D';
all_tests{4}='find shared dimension in 2D';
all_tests{5}='EEG video example';
all_tests{6}='kernel CorrCA';
all_tests{7}='kernel LDA';
all_tests{8}='motor ratings';
all_tests{10}='HBN ratings';
all_tests{11}='EEG ERP example';
all_tests{13}='regularization';
all_tests{14}='simulations';
all_tests{15}='JMLR Reviewer 1 request'; % Reviwer of JMLR
all_tests{16}='NBDT Reviewer 1 hypothetical example';
all_tests{17}='NBDT Reviewer 2 question';

savefigure = 0;

% make one figure per test. you can pick here which of these many
% tests/figures to generate. 
for figure_nr=[16]
    
    test = all_tests{figure_nr};
    figure(figure_nr);
    
    switch test
        
        case 'kernel LDA'
            
            % generate 2D data with random phase across exemplars and
            % common radii for each class. Goal is to find common radii in first component
            C=3; N=200; phirange=pi*1; phi=rand(C,N)*phirange; amp=(0:C-1)'/C*2;
            X=[]; for l=1:N, X(:,:,l)=[sin(phi(:,l)) cos(phi(:,l))].*amp; end
            
            % generate test signal on a regular grid of first 2 dims to see the mapping.
            phitest=[0:phirange/10:phirange]'; Ctest=size(phitest,1); Ntest=5; amptest=((1:Ntest)/(Ntest+1))*1.5;
            Xtest=[]; for l=1:Ntest, Xtest(:,:,l)=[sin(phitest) cos(phitest)]*amptest(l); end
            
            % Add 10% noise
            X = X+1+0.1*randn(size(X));
            Xtest=Xtest+1;
            
            % best class separation F optimized by kLDA
            [W,Fy,Y,Ytest]=klda(X,'Gaussian',2,Xtest,5);
            
            
            % show what happened
            clf; cmap=colormap('jet');
            subplot(2,2,1);
            plot(squeeze(Xtest(:,1,:))',squeeze(Xtest(:,2,:))','color',[0.7 0.7 0.7]); hold on
            plot(squeeze(Xtest(:,1,:)) ,squeeze(Xtest(:,2,:))); axis equal;
            for i=1:C, plot(squeeze(X(i,1,:)),squeeze(X(i,2,:)),'.','color',cmap(ceil(i/C*length(cmap)),:)); end
            xlabel('X_1'); ylabel('X_2'); title('Original data, X'); hold off
            subplot(2,2,2);
            plot(squeeze(Ytest(:,1,:))',squeeze(Ytest(:,2,:))','color',[0.7 0.7 0.7]); hold on
            plot(squeeze(Ytest(:,1,:)) ,squeeze(Ytest(:,2,:))); axis equal;
            for i=1:C, plot(squeeze(Y(i,1,:)),squeeze(Y(i,2,:)),'.','color',cmap(ceil(i/C*length(cmap)),:));  end
            xlabel('Y_1'); ylabel('Y_2'); title('Kernel LDA, Y=\alpha^T \phi(X)');  hold off
            subplot(2,2,3);
            plot(squeeze(Ytest(:,1,:))',repmat(phitest,1,length(amptest))','color',[0.7 0.7 0.7]); hold on
            plot(squeeze(Ytest(:,1,:)) ,repmat(phitest,1,length(amptest)) ,'.-'); hold off
            xlabel('Y_1'); ylabel('Phase')
            subplot(2,2,4);
            plot(squeeze(Ytest(:,1,:))',repmat(amptest,length(phitest),1)','color',[0.7 0.7 0.7]); hold on
            plot(squeeze(Ytest(:,1,:)) ,repmat(amptest,length(phitest),1) ,'.-'); hold off;
            xlabel('Y_1'); ylabel('Amplitude')
            
            sublabel;
            
        case 'kernel CorrCA'
            
            % generate 2D data with common phase and random radii across subjects
            % goal is to find common phase in first component
            T=40; N=2; phirange=pi*2; phi=rand(T,N)*phirange; amp=rand(T,1)+0.5; % shared phase
            X=[]; for l=1:N, X(:,:,l)=[sin(phi(:,l)) cos(phi(:,l))].*amp; end
            
            % generate test signal on a regular grid of first 2 dims to see the mapping.
            phitest=[0:phirange/20:phirange]'; Ntest=5; amptest=((1:Ntest)/(Ntest+1)+0.5);
            Xtest=[]; for l=1:Ntest, Xtest(:,:,l)=[sin(phitest) cos(phitest)]*amptest(l); end
            
            % inter subject correlation in original data
            ISCx=corrcoef([X(:,:,1),X(:,:,2)]); Dx=size(X,2); ISCx=diag(ISCx(1:Dx,Dx+1:2*Dx));
            
            % intersubject correlation optimized by kCorrCA
            [W,ISCy,Y,Ytest]=kcorrca(X,'Gaussian',1,Xtest,'mean',15);
            
            
            % show all time courses
            D=2; % just show 2
            clf, clear h;
            for i=1:D, h(i)    =subplot(2*D,3,i*3-2); plot([X(:,i,1) X(:,i,2)],'.-'); text(0.7,0.9,['ISC=' num2str(ISCx(i),1)],'units','normalized'); ylabel(['X_' num2str(i)]); end; xlabel('time (samples)')
            for i=1:D, h(i+D+1)=subplot(2*D,3,i*3-1); plot([Y(:,i,1) Y(:,i,2)],'.-'); text(0.7,0.9,['ISC=' num2str(ISCy(i),1)],'units','normalized'); ylabel(['Y_' num2str(i)]); end; xlabel('time (samples)')
            subplot(2*D,3,1); title('Original data, X');
            subplot(2*D,3,2); title('Kernel CorrCA, Y=\alpha^T \phi(X)');
            hlegend=legend('Subject 1','Subject 2'); hlegend.Position([1 2])=[0.33 0.47];
            
            % show first two dimensions that have been transformed
            h(D+1)=subplot(2,3,4);
            plot(X(:,1,1),X(:,2,1),'.',X(:,1,2),X(:,2,2),'.'); hold on
            plot(squeeze(Xtest(:,1,:)),squeeze(Xtest(:,2,:)));
            plot([X(:,1,1) X(:,1,2)]',[X(:,2,1) X(:,2,2)]','color',[0.7 0.7 0.7]); hold off
            axis equal;
            xlabel('X_1'); ylabel('X_2');  axis equal
            
            h(2*D+2)=subplot(2,3,5);
            plot(Y(:,1,1),Y(:,2,1),'.',Y(:,1,2),Y(:,2,2),'.'); hold on
            plot(squeeze(Ytest(:,1,:)),squeeze(Ytest(:,2,:)));
            plot([Y(:,1,1) Y(:,1,2)]',[Y(:,2,1) Y(:,2,2)]','color',[0.7 0.7 0.7]);  hold off
            xlabel('Y_1'); ylabel('Y_2');
            
            h(end+1)=subplot(2,3,3);
            plot(squeeze(Ytest(:,1,:))',repmat(phitest/pi*180,1,length(amptest))','color',[0.7 0.7 0.7]); hold on
            plot(squeeze(Ytest(:,1,:)) ,repmat(phitest/pi*180,1,length(amptest)) ,'.-'); hold off
            xlabel('Y_1'); ylabel('Phase (deg)')
            h(end+1)=subplot(2,3,6);
            plot(squeeze(Ytest(:,1,:))',repmat(amptest,length(phitest),1)','color',[0.7 0.7 0.7]); hold on
            plot(squeeze(Ytest(:,1,:)) ,repmat(amptest,length(phitest),1) ,'.-'); hold off;
            xlabel('Y_1'); ylabel('Amplitude')
            
            sublabel(h);
            
        case 'EEG video example'
            
            
            resultFile='./EEG-Cohen-2017/CorrCAvsMCCArocco.mat'; % MCCA not regularized
            resultFile='./EEG-Cohen-2017/CorrCAvsMCCAroccoReg12.mat'; % MCCA regularized wtih K=10
            if exist(resultFile)
                load(resultFile,'Kmcca','Kisc','sISC','srho','ISC','rho','A','Amcca','ISCtest','rhotest');
            else
                
                load('./EEG-Cohen-2017/sunday_at_roccos.mat','X');
                [T,D,N]=size(X);
                
                d=D*ones(1,N); % for mcca.m
                kReg=12*ones(1,N); % regularize
                
                % determine number of significant components for CorrCA and MCCA
                [W,ISC,Y,A]=corrca(X);
                [V,rho,Amcca]=mcca(X(:,:),d,[],kReg);
                
                % circular stats
                % want to use the same permutations for cca and mcca, so not
                % using est_K_shift here
                nPerms = 100;
                for p=1:nPerms
                    p
                    sEEG=zeros(size(X));
                    for n=1:N
                        sEEG(:,:,n)=circshift(X(:, :, n), randi(T));
                    end
                    [~,tISC,~,~]=corrca(sEEG);
                    [~,trho,~]=mcca(sEEG(:,:),d,[],kReg);
                    sISC(:,p)=tISC;
                    srho(:,p)=trho;
                end
                
                pvals = mean(repmat(rho, 1, nPerms) < repmat(srho(1, :), size(srho, 1), 1), 2);
                Kmcca = sum(pvals < 0.05);
                
                pvals = mean(repmat(ISC, 1, nPerms) < repmat(sISC(1, :), size(sISC, 1), 1), 2);
                Kisc = sum(pvals < 0.05);
                
                % now cross-validate to get the ISC comparison
                nFolds=5;
                nSamplesTrunc=floor(T/nFolds)*nFolds;
                samples2D=reshape(1:nSamplesTrunc,floor(T/nFolds),nFolds);
                
                for f=1:nFolds
                    EEGtrain=X(samples2D(:,setdiff(1:nFolds,f)),:,:);
                    EEGtest=X(samples2D(:,f),:,:);
                    
                    % train
                    [W,~,~,~]=corrca(EEGtrain);
                    % test
                    [~,tISC,Y,~]=corrca(EEGtest,'fixed',W);
                    ISCtest(:,f)=tISC;
                    
                    % train & test for M-CCA
                    [~,~,~,trhotest]=mcca(EEGtrain(:,:),d,EEGtest(:,:),kReg);
                    rhotest(:,f)=trhotest;
                    
                end
                
                % test for the best recularization constant kReg for MCCA
                if 0
                    clear muMccas
                    for kall=21:30
                        kall
                        kReg=kall*ones(1,N); % regularize
                        clear rhotest
                        for f=1:nFolds
                            EEGtrain=X(samples2D(:,setdiff(1:nFolds,f)),:,:);
                            EEGtest=X(samples2D(:,f),:,:);
                            
                            % train & test for M-CCA
                            [~,~,~,trhotest]=mcca(EEGtrain(:,:),d,EEGtest(:,:),kReg);
                            rhotest(:,f)=trhotest;
                            
                        end
                        muMccas(kall,:)=nanmean( rhotest(1:5,:),2);
                    end
                    plot(sum(muMccas,2));
                    % shows that k=12 had the highest cross validation of
                    % first 3 components (component 4 and 5 look like
                    % noise)
                end
                
                
                save(resultFile,'Kmcca','Kisc','sISC','srho','ISC','rho','A','Amcca','ISCtest','rhotest','kReg');
                
            end
            
            % draw figure
            % plot results
            clf
            nShow=5;
            hs(1)=subplot(321); hold on; box on
            hbar(1)=bar((1:nShow)-0.15,ISC(1:nShow),'barwidth',0.3);
            hbar(3)=bar((1:nShow)+0.15,rho(1:nShow),'FaceColor',[0.00 0.7 0.70],'barwidth',0.3);
            hbar(2)=bar((Kisc+1:nShow)-0.15,ISC(Kisc+1:nShow),'FaceColor',[0.75 0.75 1],'barwidth',0.3);
            hbar(4)=bar((Kmcca+1:nShow)+0.15,rho(Kmcca+1:nShow),'FaceColor',[0.75 1 1],'barwidth',0.3);
            set(gca,'Xtick',1:nShow); set(gca,'Ytick',0:0.05:0.15);
            xlabel('Component'); ylabel('ISC');
            htit(1)=title('Train');
            
            hs(2)=subplot(322);hold on; box on
            [semsIscs,muIscs]=nansem(ISCtest,2); % how come there are NaN here?
            [semsMccas,muMccas]=nansem( rhotest,2);
            hbar(1)=bar((1:nShow)-0.15,muIscs(1:nShow),'barwidth',0.3);
            hbar(3)=bar((1:nShow)+0.15,muMccas(1:nShow),'FaceColor',[0.00 0.7 0.70],'barwidth',0.3);
            hbar(2)=bar((Kisc+1:nShow)-0.15,muIscs(Kisc+1:nShow),'FaceColor',[0.75 0.75 1],'barwidth',0.3);
            hbar(4)=bar((Kmcca+1:nShow)+0.15,muMccas(Kmcca+1:nShow),'FaceColor',[0.75 1 1],'barwidth',0.3);
            herr1=errorbar((1:nShow)-0.15,muIscs(1:nShow),semsIscs(1:nShow),'r','LineStyle','none');
            herr2=errorbar((1:nShow)+0.15,muMccas(1:nShow),semsMccas(1:nShow),'r','LineStyle','none');
            hlg=legend('CorrCA','MCCA','CorrCA (n.s.)','MCCA (n.s.)');set(hlg,'Box','off');
            yl=ylim;
            ylim([0 yl(2)+0.001]);
            set(gca,'Xtick',1:nShow); set(gca,'Ytick',[0 0.025 0.05]);
            xlabel('Component');
            htit(2)=title('Test');
            
            % topoplots
            subs2show=[1 5 9 13 17];
            clims=[];
            nRows=5; rowOffset=2; nCols=numel(subs2show)+1;
            for r=rowOffset+1:nRows
                hs(r,1)=subplot(nRows,nCols,(r-1)*nCols+1);
                topoplot(A(:,r-rowOffset),'EEG-Cohen-2017//BioSemi64.loc','electrodes','off');
                clims=cat(1,clims,caxis);
                for c=2:nCols
                    hs(r,c)=subplot(nRows,nCols,(r-1)*nCols+c);
                    topoplot(Amcca{subs2show(c-1)}(:,r-rowOffset),'EEG-Cohen-2017/BioSemi64.loc','electrodes','off');
                    clims=cat(1,clims,caxis);
                end
            end
            
            % reposition scalp maps (hard-coded)
            delRight=0.03; delLeft=0.03;
            for i=rowOffset+1:nRows
                pos=get(hs(i,2),'Position');
                set(hs(i,2),'Position',[pos(1)+delRight pos(2) pos(3) pos(4)]);
                
                pos=get(hs(i,3),'Position');
                set(hs(i,3),'Position',[pos(1)+delRight/2 pos(2) pos(3) pos(4)]);
                
                pos=get(hs(i,5),'Position');
                set(hs(i,5),'Position',[pos(1)-delLeft/2 pos(2) pos(3) pos(4)]);
                
                pos=get(hs(i,6),'Position');
                set(hs(i,6),'Position',[pos(1)-delLeft pos(2) pos(3) pos(4)]);
            end
            
            % add titles
            set(get(hs(3,1),'Title'),'String','CorrCA');
            for c=2:nCols
                set(get(hs(rowOffset+1,c),'Title'),'String',['MCCA: S' num2str(subs2show(c-1))]);
            end
            
            % add dreaded colorbar
            hcb=colorbar('east');
            %set(hs(rowOffset+1,:),'Clim',[min(clims(:,1)) max(clims(:,2))]); %looks ugly
            cbPos=get(hcb,'Position');
            set(hcb,'Position',[0.25 cbPos(2) cbPos(3) cbPos(4)]);
            clim=get(hcb,'limits'); set(hcb,'Ticks',[clim(1) 0 clim(2)]);
            set(hcb,'TickLabels',{'min','0','max'});
            
            sublabel([hs(1,1) hs(1,2) hs(3,1) hs(3,2)],-16,-20,'FontSize',16);
            
        case 'EEG ERP example'
            
            % ERN data from Parra, Lucas, et al. "Linear spatial
            % integration for single-trial detection in encephalography."
            % NeuroImage 17.1 (2002): 223-230.
            %
            % back in 2002, this is the messy way I managed data. To my
            % defence, I got the data this way.
            datafile = 'ERN-Parra-2002/4_02_errors(46).bin'; N=46; D=64; fs=250; T=176;
            locfile  = 'ERN-Parra-2002/coord.loc';
            
            clf
            
            time=1000*((1:T)/fs-0.20);
            fid = fopen(datafile,'r','b');
            error = fread(fid,[128,T*N],'int16');
            error = error(1:D,:);
            error = reshape(error,[D,T,N]);
            fclose(fid);
            
            X=permute(error,[2 1 3]);
            
            % compute inter-trial correlation on the raw data
            [~,ITCx]=corrca(X,'fixed',eye(D));
            
            subplot(5,3,1);
            SNRx=(ITCx+1/(N-1))./(1-ITCx)*sqrt(N);
            bar(SNRx); xlabel('EEG channels'); ylabel('SNR*\surd N')
            axis([0 D 0 3.5])
            
            % electrodes picked 'by hand' by looking at ERP (FCz and Cz)
            [~,indx]=sort(SNRx);
            n=4; for i=[7 22]
                x=squeeze(X(:,i,:));
                subplot(5,3,n); n=n+3; ploterp(x',[],0,time);
                title(['trial averaged x_{' num2str(i) '}(t)']);
                ylabel('\mu V')
                subplot(5,3,n); n=n+3; imagesc(time,[],x');
                ylabel('trials');
                title(['single trial x_{' num2str(i) '}(t)'])
                caxis([-1 1]*prctile(abs(x(:)),100))
            end
            xlabel('time (ms)');
            
            % find corrca components, and evalute performance on
            % leave-one-out sample
            for i=N:-1:1
                Xtrain = X; Xtrain(:,:,i)=[];
                [W(:,:,i),~,~,A(:,:,i)] = corrca(Xtrain,'shrinkage',0.4);
                Sign = -diag(sign(diag(inv(W(:,:,N))*W(:,:,i)))); % fix arbitrary sign
                Ytest(:,:,i) = X(:,:,i)*W(:,:,i)*Sign;
            end
            [~,ITCy,~,~]=corrca(Ytest,'fixed',eye(D));
            
            
            % find corrca components + establish significance through
            % random circular shifting
            [K, p] = est_K_shift(X, 100, 0.4);
            disp(['Number of significant components using circular shuffle p<0.05:' num2str(K)]);
            
            subplot(5,3,2);
            SNRy=(ITCy+1/(N-1))./(1-ITCx)*sqrt(N);
            bar(SNRy); hold on
            i=K+1:D; bar(i,SNRy(i),'FaceColor',[0.75 0.75 1]);
            axis([0 D 0 3.5])
            xlabel('CorrCA components'); ylabel('SNR*\surd N ');
            
            n=5; for i=1:2
                y=squeeze(Ytest(:,i,:));
                subplot(5,3,n); n=n+3;
                ploterp(y',[],0,time);
                title(['trial averaged y_' num2str(i) '(t)'])
                ylabel('\mu V')
                subplot(2,3,3*i);
                topoplot(A(:,i),locfile,'electrodes','numbers');
                title(['a_' num2str(i)])
                h=subplot(5,3,n); n=n+3;
                imagesc(time,[],y');
                ylabel('trials');
                title(['single trial y_' num2str(i) '(t)'])
                caxis([-1 1]*prctile(abs(y(:)),100))
                pos=get(h,'Position');
            end
            xlabel('time (ms)');
            h=subplot(5,3,14);
            drawnow
            pos=get(h,'Position');
            hcbar=colorbar; title(hcbar,'\mu V')
            set(h,'Position',pos);
            
            % add some helpful lines
            for i=[4 5 7 8]
                subplot(5,3,i);
                ax = axis;
                hold on;
                plot([0 0], ax(3:4), 'k')
                plot([100 100], ax(3:4), 'r')
            end
            
            sublabel([],-10,-30);
            
        case 'CorrCA version speed'
            
            T=5000;D=64;N=15; % samples/classes, dimensions, subjects/examplars
            X = randn(T,D,N)+repmat(randn(T,D),1,1,N); % add same mean time course to all subjects
            
            tic; [Wcca1,ISC1]=corrca(X,'version',1); t=toc; disp(['CorrCA vs 1: ' num2str(t,3) ' seconds'])
            tic; [Wcca2,ISC2]=corrca(X,'version',2); t=toc; disp(['CorrCA vs 2: ' num2str(t,3) ' seconds'])
            
            % show that the result is the same
            subplot(2,2,1); imagesc(abs(inv(Wcca1)*Wcca2)); axis square; caxis([0 1])
            title('W_1^{-1} W_2')
            subplot(2,2,2); plot([ISC1./ISC2]); axis tight
            title('ISC_1/ISC_2')
            
            sublabel
            
        case 'find shared dimension in 2D'
            
            % 2D matrix to rotate by alpha and scale
            alpha1=pi/4;alpha2=pi/3;
            M = [sin(alpha1) cos(alpha1);cos(alpha2) -sin(alpha2)]*diag([1 1]);
            
            
            % generate data with 1 dimension shared across 2 subjects
            x=randn(20,3); % sharing first dimension + 2 dimensions each per subject
            X=cat(3,x(:,[1 2])*M,x(:,[1 3])*M); % rotating them the same accross subjects
            [T,D,N]=size(X); % samples, dimensions, subjects
            
            % what if the noise direction is not the same?
            %             alpha0=0; alpha1=pi/4;alpha2=pi/3;
            %             M1=[sin(alpha0) cos(alpha0);cos(alpha1) -sin(alpha1)];
            %             M2=[sin(alpha0) cos(alpha0);cos(alpha2) -sin(alpha2)];
            %             X=cat(3,x(:,[1 2])*M1,x(:,[1 3])*M2);
            % it still finds the common direction, but there is noise to
            % it. Basically acts the way regular additive noise would
            % affect results.
            
            % Add 5% noise
            X = X+0.0*randn(size(X));
            
            % inter subject correlation in original data
            ISCx=corrcoef([X(:,:,1),X(:,:,2)]); ISCx=diag(ISCx(1:D,D+1:2*D));;
            
            % intersubject correlation optimized by CorrCA
            %             [W,ISCy,Y]=corrca(X,'version',2);
            [K, p, W, ISCy, Y, A] = est_K_shift(X, 1000, 0);
            
            % show all time courses
            clf, clear h;
            for i=1:D, h(i)    =subplot(2*D,2,i*2-1); plot(X(:,i,1),'v-'); hold on; plot(X(:,i,2),'^-'); text(0.6,0.87,['ISC=' num2str(ISCx(i),1)],'units','normalized','backgroundcolor',[1 1 1]); ylabel(['X_' num2str(i)]); end; xlabel('time (samples)')
            for i=1:D, h(i+D+1)=subplot(2*D,2,i*2  ); plot(Y(:,i,1),'v-'); hold on; plot(Y(:,i,2),'^-'); text(0.6,0.87,['ISC=' num2str(ISCy(i),1) ', p=' num2str(p(i),1)],'units','normalized','backgroundcolor',[1 1 1]); ylabel(['Y_' num2str(i)]); end; xlabel('time (samples)')
            subplot(2*D,2,1); title('Original data, X');
            subplot(2*D,2,2); title('Correlated Components, Y=W'' X');
            hlegend=legend('Subject 1','Subject 2'); hlegend.Position([1 2])=[0.45 0.476];
            
            % show first two dimensions that have been rotated as scatter plot
            h(D+1)=subplot(2,2,3);
            
            plot(X(:,1,1),X(:,2,1),'v',X(:,1,2),X(:,2,2),'^'); hold on
            plot([X(:,1,1) X(:,1,2)]',[X(:,2,1) X(:,2,2)]','color',[0.8 0.8 0.8]); hold off
            axis equal;
            xlabel('X_1'); ylabel('X_2')
            
            if 1 % show the MCCA results
                [V,rho,Amcca]=mcca(X(:,:),[2 2]);
                V = {V(1:2,1:2),V(3:4,1:2)};
                arrow([0 0], Amcca{1}(:,1), 'color',[1 0.5 0.5]);
                arrow([0 0], Amcca{1}(:,2), 'color',[0.5 0.5 1]);
                arrow([0 0], Amcca{2}(:,1), 'color',[1 0.5 0.5 ]);
                arrow([0 0], Amcca{2}(:,2), 'color',[0.5 0.5 1 ]);
                disp(['Deviation of component 1 between CorrCA and MCCA:'  ...
                    num2str([
                    subspace(A(:,1), Amcca{1}(:,1)) , ...
                    subspace(A(:,1), Amcca{2}(:,1))]/pi*180,2)])
            end
            %
            arrow([0 0], A(:,1), 'color','r');
            arrow([0 0], A(:,2), 'color','b');
            
            h(2*D+2)=subplot(2,2,4);
            plot(Y(:,1,1),Y(:,2,1),'v',Y(:,1,2),Y(:,2,2),'^'); hold on
            plot([Y(:,1,1) Y(:,1,2)]',[Y(:,2,1) Y(:,2,2)]','color',[0.8 0.8 0.8]);
            axis equal; hold off
            xlabel('Y_1'); ylabel('Y_2')
            Ydirection =A'*W;
            arrow([0 0], Ydirection(:,1), 'color','r')
            arrow([0 0], Ydirection(:,2), 'color','b')
            
            sublabel(h);
            
            
        case 'find shared dimension in 3D'
            
            % 3D matrix to rotate first 2 dimension by alpha and scale all 3
            alpha=pi/4;
            M = [sin(alpha) cos(alpha) 0;cos(alpha) -sin(alpha) 0; 0 0 1]*diag([2 1 0.5]);
            
            % generate data with 1 dimension shared across 3 subjects
            x=randn(20,7); % sharing first dimension + 2 dimensions each per subject
            X=cat(3,x(:,[1 2 3])*M,x(:,[1 4 5])*M,x(:,[1 6 7])*M); % rotating them the same accross subjects
            [T,D,N]=size(X); % samples, dimensions, subjects
            
            % Add 5% noise
            X = X+0.0*randn(size(X));
            
            % inter subject correlation in original data
            ISCx=corrcoef([X(:,:,1),X(:,:,2)]); ISCx=diag(ISCx(1:D,D+1:2*D));;
            
            % intersubject correlation optimized by CorrCA
            [W,ISCy,Y]=corrca(X,'version',2);
            
            % show all time courses
            clf
            for i=1:D, h(i)  =subplot(2*D,2,i*2-1); plot([X(:,i,1) X(:,i,2)],'*-'); text(0.7,0.9,['ISC=' num2str(ISCx(i),1)],'units','normalized'); ylabel(['X_' num2str(i)]); end; xlabel('time (samples)')
            for i=1:D, h(i+D+1)=subplot(2*D,2,i*2  ); plot([Y(:,i,1) Y(:,i,2)],'*-'); text(0.7,0.9,['ISC=' num2str(ISCy(i),1)],'units','normalized'); ylabel(['Y_' num2str(i)]); end; xlabel('time (samples)')
            subplot(2*D,2,1); title('Original data, X');
            subplot(2*D,2,2); title('Correlated Components, Y=W^T X');
            
            % show first two dimensions that have been rotated as scatter plot
            h(D+1)=subplot(2,2,3);
            plot(X(:,1,1),X(:,2,1),'*',X(:,1,2),X(:,2,2),'*'); hold on
            plot([X(:,1,1) X(:,1,2)]',[X(:,2,1) X(:,2,2)]','color',[0.8 0.8 0.8]); hold off
            axis equal; legend('Subject 1','Subject 2')
            xlabel('X_1'); ylabel('X_2')
            
            h(2*D+2)=subplot(2,2,4);
            plot(Y(:,1,1),Y(:,2,1),'*',Y(:,1,2),Y(:,2,2),'*'); hold on
            plot([Y(:,1,1) Y(:,1,2)]',[Y(:,2,1) Y(:,2,2)]','color',[0.8 0.8 0.8]); axis equal; hold off
            legend('Subject 1','Subject 2');
            xlabel('Y_1'); ylabel('Y_2')
            
            sublabel(h);
            
            
        case 'CorrCA vs FLD vs JD' % show that the two give the same result, provided zero mean per subject
            
            clf
            T=10;D=5;N=2; % samples/classes, dimensions, subjects/examplars
            X = randn(T,D,N)+repmat(randn(T,D),1,1,N); % add same mean time course to all subjects
            
            [Wcca,ISC]=corrca(X);
            [Wlda,F]=lda(X);
            [Wjd,SNR]=corrca(X,'version',3);
            
            % show that the result differs
            subplot(2,2,1); imagesc(abs(inv(Wlda)*Wcca)); axis square; caxis([0 1])
            title('Non-zero mean: W_{LDA}^{-1} W_{CorrCA}')
            subplot(2,2,3); imagesc(abs(inv(Wjd)*Wcca)); axis square; caxis([0 1])
            title('Non-zero mean: W_{JD}^{-1} W_{CorrCA}')
            
            % subtract mean per subject
            for l=1:N, X(:,:,l)=X(:,:,l)-repmat(mean(X(:,:,l)),T,1); end
            
            [Wcca,ISC]=corrca(X);
            [Wlda,F]=lda(X);
            [Wjd,SNR]=corrca(X,'version',3);
            
            % show that the result is the same
            subplot(2,2,2); imagesc(abs(inv(Wlda)*Wcca)); axis square; caxis([0 1])
            title('Zero mean: W_{LDA}^{-1} W_{CorrCA}')
            subplot(2,2,4); imagesc(abs(inv(Wjd)*Wcca)); axis square; caxis([0 1])
            title('Zero mean: W_{JD}^{-1} W_{CorrCA}')
            
            sublabel;
            
        case 'motor ratings'
            
            % Data courtesy of Valentin Rousson: Rousson, Valentin, Gasser,
            % Theo, Caflisch, Jon and Largo, Remo(2006) 'Reliability of the
            % Zurich Neuromotor Assessment',The Clinical Neuropsychologist,
            % 22:1,60-72
            
            clear X Xtest
            datapath='Rousson_child_motor_assessment/';
            [child, age, sex, task, rater1, rater2, rater2b, rater1_sds, rater2_sds, rater2b_sds] = ...
                textread([datapath 'inter.rater.dat'],'%d%d%d%d%f%f%f%f%f%f\n','headerlines',1);
            labels = importdata([datapath 'inter.rater.labels']);
            for i=unique(child)'
                for j=unique(task)'
                    indx=find(child==i & task==j);
                    X(i,j,:)     = [rater1_sds(indx); rater2_sds(indx)];
                    Xtest(i,j,:) = [rater2_sds(indx); rater2b_sds(indx)];
                end
            end
            tasks=1:12;
            X=X(:,tasks,:);Xtest=Xtest(:,tasks,:);
            labels = {labels{tasks}};
            
            [T,D,N]=size(X); % samples, dimensions, raters
            
            % inter-rater correlation in original data
            [~,IRCx,    ~,~,px    ]=corrca(X,    'fixed',eye(D));
            [~,IRCxtest,~,~,pxtest]=corrca(Xtest,'fixed',eye(D));
            
            % inter-rater correlation optimized with CorrCA
            [W,IRCy,Y,A]=corrca(X);
            
            % inter-rater correlation in component space for test data
            [~,IRCytest,Ytest,~,pytest]=corrca(Xtest,'fixed',W);
            
            % number of significant components using F-stats with
            % Bonferroni correction on test data
            Ktest = sum(pytest*D<0.05);
            
            % p values and number of significant components on training data
            [Ktrain, p] = est_K_shift(X, 1000, 0);
            
            
            tmp=permute(X-repmat(mean(X),T,1),[2 1 3]); ccx=corrcoef(tmp(:,:)');
            tmp=permute(Y-repmat(mean(Y),T,1),[2 1 3]); ccy=corrcoef(tmp(:,:)');
            tmp=permute(Ytest-repmat(mean(Ytest),T,1),[2 1 3]); ccytest=corrcoef(tmp(:,:)');
            
            clf;
            tmp=[ccx(:);ccy(:);ccytest(:)]; clim=[min(tmp) max(tmp)];
            subplot(2,4,1); imagesc(ccx); caxis(clim); xlabel('Tasks'); ylabel('Tasks'); title('R_W^x')
            h=colorbar; h.Position(1)=0.05; ylabel(h,'corr. coef')
            subplot(2,4,2); bar(IRCx); hold on
            i=find(px*D>0.05); bar(i,IRCx(i),'FaceColor', [0.75 0.75 1])
            axis tight; ylim([0 1]); xlabel('Tasks'); ylabel('IRC');     title('Inter-rater')
            subplot(2,4,3); bar(IRCxtest); hold on
            i=find(pxtest*D>0.05); bar(i,IRCxtest(i),'FaceColor', [0.75 0.75 1])
            axis tight; ylim([0 1]); xlabel('Tasks'); ylabel('IRC'); title('Test-retest')
            subplot(2,4,5); imagesc(ccy); caxis([0 1]); xlabel('Components');  ylabel('Components');  title('R_W^y')
            h=colorbar; h.Position(1)=0.05; ylabel(h,'corr. coef')
            subplot(2,4,6);
            bar(IRCy); hold on
            i=Ktrain+1:D; bar(i,IRCy(i),'FaceColor', [0.75 0.75 1])
            axis tight; ylim([0 1]); xlabel('Components');  ylabel('IRC');
            
            subplot(2,4,7);
            bar(IRCytest); hold on;
            i=Ktest+1:D; bar(i,IRCytest(i),'FaceColor', [0.75 0.75 1])
            axis tight; ylim([0 1]); xlabel('Components');  ylabel('IRC');
            
            subplot(1,4,4); imagesc(W(:,1:Ktrain)); xlabel('Components'); title('V'); caxis([-1 1]*max(abs(W(:))));
            % subplot(1,4,4); imagesc(A(:,1:Ktrain)); xlabel('Components'); title('A'); caxis([-1 1]*max(abs(W(:))));
            set(gca,'ytick',[1:D],'yticklabels',labels,'YAxisLocation', 'right','YTickLabelRotation',-45)
            ylabel('Tasks'); h=get(gca,'Ylabel'); h.Position(1)=-0.5;
            set(gca,'YTickLabelRotation',-45)
            h=colorbar('southoutside');
            pos=get(gca,'position');pos([2 4])=[0.19 0.73]; set(gca,'position',pos);
            h.Position(2)=0.05;
            
            sublabel([],-10,-30);
            
            
        case 'HBN ratings'
            
            % Data courtesy of Mike Milham. Please refer to:
            % https://www.biorxiv.org/content/early/2017/07/17/149369
            
            datapath = 'HealthyBrainNetwork_Phenotypic_Data/';
            clear ratings indx
            [~,data{1}]=xlsread([datapath 'Alabama Parenting Questionnaire.xlsx'],'B2:AX735');
            indx(1,:)=[1:42 43:49];
            [~,data{2}]=xlsread([datapath 'Alabama Parenting Questionnaire.xlsx'],'BL2:DR735');
            indx(2,:)=[1 3:5 7:9 11:12 14:15 17:19 21 23:27 29:34 36:51 52:58];
            [~,demographic]=xlsread([datapath 'Alabama Parenting Questionnaire.xlsx'],'EF2:EG735');
            raw_scores = 1:42; agregate_scores=43:49;
            D = size(indx,2); T=size(data{1},1);
            ratings=zeros(T,D,2); demos=zeros(T,size(demographic,2));
            for i=1:T
                for j=1:size(demographic,2)
                    str = demographic{i,j};
                    if ~isempty(str) || ~strcmp(str,'.'), demos(i,j)=str2num(str); else demos(i,j) = NaN; end
                end
                for r=1:2
                    for j=1:D
                        str = data{r}{i,indx(r,j)};
                        if isempty(str)  % missing entry
                            ratings(i,j,r)=NaN;
                        else
                            if strcmp(str,'.')  % data not collected yet
                                ratings(i,j,r)=NaN;
                            else
                                ratings(i,j,r) = str2num(str);
                            end
                        end
                    end
                end
            end
            % removed missing data
            indx = find(isnan(sum(sum(ratings,3),2)));
            ratings(indx,:,:)=[];
            demos(indx,:)=[];
            
            % inter-rater correlation of original agregate ratings
            Xall=ratings(:,raw_scores,:);
            [Tall,D,N]= size(Xall);
            [~,IRCraw,~,~,pxall]=corrca(Xall,'fixed',eye(D));
            [~,IRCagregate]=corrca(ratings(:,agregate_scores,:),'fixed',eye(length(agregate_scores)));
            [~,p,~,stats]= ttest2(IRCraw,IRCagregate);
            disp(['raw scores     =' num2str(mean(IRCraw),2)      '+/-' num2str(std(IRCraw)     ,2)])
            disp(['agregate scores=' num2str(mean(IRCagregate),2) '+/-' num2str(std(IRCagregate),2)])
            disp(['Do raw and clinical agregate scores differ? 2-sample t-test, t(' num2str(stats.df) ')=' num2str(stats.tstat,2) ', p=' num2str(p,2)])
            
            % correlation of questions in original data
            [~,IRCxall]=corrca(Xall,'fixed',eye(D));
            tmp=permute(Xall-repmat(mean(Xall),Tall,1),[2 1 3]); ccxall=corrcoef(tmp(:,:)');
            Nrand=100; clear Wall IRCytrain IRCytest IRCytrainK IRCytestK
            for nrand = 1:Nrand
                
                % split data into train and test at random
                indx = randperm(Tall); T2=round(Tall/2);
                itrain = indx(1:T2);
                itest  = indx(T2+1:Tall);
                % now sample with replacement from each half (to get
                % bootstrap estimates of STD)
                itrain = itrain(randi(T2,T2,1));
                itest = itest(randi(T2,T2,1));
                Xtrain = ratings(itrain,raw_scores,:);
                Xtest  = ratings(itest ,raw_scores,:);
                Ttrain= size(Xtrain,1);
                Ttest = size(Xtest,1);
                
                % inter-rater correlation optimized with CorrCA
                [Wall(:,:,nrand),IRCytrain(:,nrand),Ytrain,Aall(:,:,nrand)]=corrca(Xtrain,'shrinkage',0.5);
                
                % inter-rater correlation in component space for test data
                [~,IRCytest(:,nrand),Ytest,~,pytest]=corrca(Xtest,'fixed',Wall(:,:,nrand));
                Ktest_F(nrand) = sum(pytest*D < 0.05);
                
                % try with kCorrCA
                [~,IRCytrainK(:,nrand),~,~,IRCytestK(:,nrand)]=kcorrca(Xtrain,'Gaussian',50,Xtest,'mean',10);
                
                % try with mcca
                [~,IRCytrainM(:,nrand),~,IRCytestM(:,nrand)]=mcca(Xtrain(:,:),[D D],Xtest(:,:));
                
            end
            
            % Report number of sinificant componets on test sets
            Ktest_F = floor(median(Ktest_F));
            disp(['Significant components on test data (meadian over ' ...
                num2str(Nrand) ' train/test splits) using F-statistic: ' ...
                num2str(Ktest_F) ]);
            
            % Now also establish significance through random circular
            % shifting using all the data as training data
            [Ktrain_shift, p] = est_K_shift(ratings(:,raw_scores,:), 1000, 0.5);
            disp(['Number of significant components using circ shift test on all data (p<0.05): ' num2str(Ktrain_shift)])
            
            % do the same with the kernel variant (Stefan, can you put this
            % in a separate function as in est_K_shift.m ?
            Nshuffle = 100;
            [~, r_] = kcorrca(ratings(:,raw_scores,:), 'Gaussian', 50, [], 'mean', 10);
            clear r_surro
            data_surro = ratings(:,raw_scores,:);
            for ishuffle = 1:Nshuffle
                for in = 1%:N
                    data_surro(:, :, in) = circshift(ratings(:, raw_scores, in), randi(Tall));
                end
                [~, r__] = kcorrca(data_surro, 'Gaussian', 50, [], 'mean', 10);
                r_surro(:, ishuffle) = r__;
            end
            pkernel = mean(repmat(r_, 1, Nshuffle) < repmat(r_surro(1, :), size(r__, 1), 1), 2);
            Kkernel = sum(pkernel < 0.05);
            disp(['Number of significant components for kernel-CorrCA, using circ shift test on all data (p<0.05): ' ...
                num2str(Kkernel)])
            
            %  MCCA circular shift statistics
            Nshuffle = 100;
            tmp = ratings(:,raw_scores,:);
            [~, r_] = mcca(tmp(:,:),[D D]);
            clear r_surro
            data_surro = ratings(:,raw_scores,:);
            for ishuffle = 1:Nshuffle
                for in = 1:N
                    data_surro(:, :, in) = circshift(ratings(:, raw_scores, in), randi(Tall));
                end
                [~,r__]=mcca(data_surro(:,:),[D D]);
                r_surro(:, ishuffle) = r__;
            end
            pmcca = mean(repmat(r_, 1, Nshuffle) < repmat(r_surro(1, :), size(r__, 1), 1), 2);
            Kmcca = sum(pmcca < 0.05);
            disp(['Number of significant components for MCCA, using circ shift test on all data (p<0.05): ' ...
                num2str(Kmcca)])
            
            disp(['ISC summed over K=' num2str(Kmcca) ' components for CorrCA: ' num2str(sum(mean(IRCytest (1:4,:)')),2)])
            disp(['ISC summed over K=' num2str(Kmcca) ' components for   MCCA: ' num2str(sum(mean(IRCytestM(1:4,:)')),2)])
            
            tmp=permute(Ytrain-repmat(mean(Ytrain),Ttrain,1),[2 1 3]); ccytrain=corrcoef(tmp(:,:)');
            tmp=permute(Ytest -repmat(mean(Ytest ),Ttest ,1),[2 1 3]); ccytest =corrcoef(tmp(:,:)');
            
            % now show some results
            clf;
            subplot(2,4,1);
            i=find(pxall*D<0.05);  bar(i,IRCxall(i)); hold on
            i=find(pxall*D>=0.05); bar(i,IRCxall(i),'FaceColor',[0.75 0.75 1]);
            axis tight;
            ylim([-0 1]); xlabel('Questions'); ylabel('IRC'); title('All X')
            
            subplot(2,4,2);
            bar(mean(IRCytrain,2)); hold on
            i = Ktrain_shift+1:D; bar(i,mean(IRCytrain(i,:),2),'FaceColor',[0.75 0.75 1]);
            errorbar(mean(IRCytrain'),std(IRCytrain'),'.r'); hold off
            axis tight; ylim([-0 1]); xlim([0 D+1]); xlabel('Components');  ylabel('IRC'); title('Train Y')
            legend('CorrCA')
            
            subplot(2,4,3);
            bar(mean(IRCytest,2)); hold on
            i = Ktest_F+1:D; bar(i,mean(IRCytest(i,:),2),'FaceColor',[0.75 0.75 1]);
            errorbar(mean(IRCytest'),std(IRCytest'),'.r'); hold off
            axis tight; ylim([-0 1]); xlabel('Components');  ylabel('IRC'); title('Test Y')
            legend('CorrCA')
            
            clim=[min(ccxall(:)) max(ccxall(:))];
            subplot(2,4,5); imagesc(ccxall); caxis(clim); xlabel('Questions'); ylabel('Questions');
            h=colorbar; h.Position(1)=0.05; ylabel(h,'corr. coef')
            
            subplot(2,4,6);
            KD=size(IRCytrainK,1);
            i = 1:Kmcca; bar(i-0.15,mean(IRCytrainM(i,:),2),'FaceColor',[0.00 0.7 0.70],'barwidth',0.3); hold on
            i = 1:Kkernel; bar(i+0.15,mean(IRCytrainK(i,:),2),'FaceColor',[1.00 1.00 0],'barwidth',0.3);
            i = Kmcca+1:KD; bar(i-0.15,mean(IRCytrainM(i,:),2),'FaceColor',[0.75 1 1],'barwidth',0.3);
            i = Kkernel+1:KD; bar(i+0.15,mean(IRCytrainK(i,:),2),'FaceColor',[1.00 1.00 0.75],'barwidth',0.3);
            errorbar([(1:KD)'-0.15,      (1:KD)'+0.15], ...
                [mean(IRCytrainM(1:KD,:)'); mean(IRCytrainK(1:KD,:)')]', ...
                [ std(IRCytrainM(1:KD,:)');  std(IRCytrainK(1:KD,:)')]','.r'); hold off
            axis tight; ylim([-0 1]); xlabel('Components');  ylabel('IRC');
            legend('MCCA','KCorrCA')
            
            subplot(2,4,7);
            bar([mean(IRCytest(1:4,:)'); mean(IRCytestM(1:4,:)'); mean(IRCytestK(1:4,:)')]'); hold on;
            errorbar([(1:4)'-0.225,       (1:4)',                  (1:4)'+0.225], ...
                [mean(IRCytest(1:4,:)'); mean(IRCytestM(1:4,:)'); mean(IRCytestK(1:4,:)')]', ...
                [ std(IRCytest(1:4,:)');  std(IRCytestM(1:4,:)');  std(IRCytestK(1:4,:)')]','.r'); hold off
            axis tight; ylim([-0 1]); xlabel('Components');  ylabel('IRC');
            legend('CorrCA','MCCA','KCorrCA')
            
            for i=1:Nrand, tmp = inv(Wall(:,:,1))*Wall(:,:,i); sgn(i)=sign(tmp(1)); end
            Wall = squeeze(Wall(:,1,:))*diag(sgn);
            %             subplot(1,4,4); barh(mean(Wall,2)); hold on
            %             errorbar(mean(Wall,2),1:D,std(Wall,[],2),'.r','horizontal'); hold off
            %             title('V'); ylim([0 D+1]); grid on
            [~,indx] = sort(abs(mean(Wall,2)),'descend');
            disp('Backward model weights for top 5 questions in component 1:')
            disp(num2str([mean(Wall(indx(1:6),:),2) indx(1:6) ],2))
            
            Aall = squeeze(Aall(:,1,:))*diag(sgn);
            subplot(1,4,4); barh(mean(Aall,2)); hold on
            errorbar(mean(Aall,2),1:D,std(Aall,[],2),'.r','horizontal'); hold off
            title('A'); ylim([0 D+1]); grid on
            [~,Aindx] = sort(abs(mean(Aall,2)),'descend');
            disp('Forward  model weights for top 5 questions in component 1:')
            disp(num2str([mean(Aall(Aindx(1:6),:),2) Aindx(1:6) ],2))
            
            
            
            
            age=demos(:,1); sex=demos(:,2);
            [r,p]=corrcoef(age,ratings(:,raw_scores,1)*mean(Wall,2));
            disp(['Correlation of component 1 with age (parent): ' num2str(r(2),2) ', p=' num2str(p(2),2)]);
            [r,p]=corrcoef(age,ratings(:,raw_scores,2)*mean(Wall,2));
            disp(['Correlation of component 1 with age (child): ' num2str(r(2),2) ', p=' num2str(p(2),2)]);
            
            set(gca,'YAxisLocation','right'); ylim([0 D+1])
            ylabel('Questions'); % h=get(gca,'Ylabel'); h.Position(1)=-0.5;
            
            sublabel([],-10,-30);
            
            
            
        case 'regularization'
            
            dataFilename='./EEG-Dmochowski-2014/superbowlData.mat';
            locfile='./EEG-Dmochowski-2014/BioSemi64.loc';
            load(dataFilename,'dataTrain','dataTest');
            nComp2Test=3;
            
            [T,D,N]=size(dataTrain);
            
            % compute performance of SVD regularization
            Kspace=nComp2Test:D;
            totalISC =zeros(numel(Kspace),2);
            Aall=zeros(D,nComp2Test,numel(Kspace));
            for k=1:numel(Kspace)
                [W,ISCtrain,Y,A]=corrca(dataTrain,'tsvd',Kspace(k));
                Aall(:,:,k)=A(:,1:nComp2Test);
                [~,ISCtest]=corrca(dataTest,'fixed',W);
                totalISC(k,:)=[sum(ISCtrain(1:nComp2Test)) sum(ISCtest(1:nComp2Test))];
            end
            
            % compute performance of shrinkage regularization
            gammaSpace=linspace(0,1,51);
            totalISCShr =zeros(numel(gammaSpace),2);
            Aallshr=zeros(D,nComp2Test,numel(gammaSpace));
            for g=1:numel(gammaSpace)
                [W,ISCtrain,Y,A]=corrca(dataTrain,'shrinkage',gammaSpace(g));
                Aallshr(:,:,g)=A(:,1:nComp2Test);
                [~,ISCtest]=corrca(dataTest,'fixed',W);
                totalISCshr(g,:)=[sum(ISCtrain(1:nComp2Test)) sum(ISCtest(1:nComp2Test))];
            end
            
            % show results
            figure(figure_nr); clf
            
            kvals2show=[62 18];
            gamma2show=[21];
            signs=ones(nComp2Test,numel(kvals2show));
            signs(1,2)=-1; signs(1,1)=1; signs(2,1)=-1; signs(2,2)=1;
            signsshr=ones(nComp2Test,numel(kvals2show));
            signsshr(1,2)=-1; signsshr(2,2)=-1; signsshr(3,2)=-1;
            nRows=nComp2Test;
            nCols=(numel(kvals2show)+numel(gamma2show))*2; % leave half for ISC plot
            for row=1:nRows
                for col=1:nCols/3
                    % no reg + tsvd
                    hs(row,col)=subplot(nRows,nCols,(row-1)*nCols+col);
                    topoplot(signs(row,col)*squeeze(Aall(:,row,kvals2show(col))),locfile,'numcontour',0,'electrodes','off');
                    if row==1
                        if col==1
                            title('No Reg.','FontWeight','normal');
                        else
                            title(['TSVD, K=' num2str(Kspace(kvals2show(col)))],'FontWeight','normal');
                        end
                    end
                    
                    % shrinkage
                    hs(row,numel(kvals2show)+1)=subplot(nRows,nCols,(row-1)*nCols+numel(kvals2show)+1);
                    topoplot(signsshr(row,col)*squeeze(Aallshr(:,row,gamma2show(1))),locfile,'numcontour',0,'electrodes','off');
                    if row==1
                        title(['Shrink, \gamma=' num2str(gammaSpace(gamma2show))], 'FontWeight','normal');
                    end
                end
            end
            h(1)=hs(1,1); % this is where we place label
            
            h(2)=subplot(222);
            plot(Kspace,totalISC);
            xlabel('TSVD, K'); ylabel('sum ISC'); box off; axis tight
            set(gca,'xdir','reverse')
            legend('train','test','location','northeast'); legend(gca,'boxoff')
            grid on
            h(3)=subplot(224);
            plot(gammaSpace,totalISCshr);
            xlabel('Shrinkage, \gamma'); ylabel('sum ISC'); box off; axis tight
            legend('train','test','location','northeast');  legend(gca,'boxoff')
            grid on
            
            sublabel(h,-10,-30);
            
        case 'simulations'
            
            cd('stefan/sims')
            addpath('..'); % for simluation plots
            addpath('../libs/export_fig'); % for figure printing
            
            plot_sim1a
            plot_sim1b
            plot_sim2a
            plot_sim2b
            plot_sim3a
            plot_sim3b
            plot_sim5a
            plot_sim5b
            plot_sim6
            plot_sim8a
            plot_sim8b
            cd('../..')
            
        case 'JMLR Reviewer 1 request'
            
            addpath('stefan')
            addpath('stefan/libs/export_fig')
            
            % Nrep: number of repetitions to obtain confidence interval
            Nrep = 100;
            
            % subjects, dimension, samples
            N = 5; D = 30; T = 200;
            
            % SNR
            db_snr = 0;
            snr = 10.^(db_snr/20)./(1+10.^(db_snr/20)); %=0.5
            
            % K: number of correlated components
            Ks = 1:3:D;
            
            clear ISC ISC_av_pca ISC_av_test_pca
            for ik = 1:length(Ks)
                for irep = 1:Nrep
                    
                    K = Ks(ik);
                    
                    [data, Y_true, A_true] = gen_data(2*T, D, N, ones(1, K), snr, 0, 1, 1, 0, 'Gauss');
                    
                    data_test = data(T+1:end, :, :);
                    data = data(1:T, :, :);
                    
                    [W, ISC(:,ik,irep)] = pca_mean(data);
                    [~, ISC_test(:,ik,irep)] = pca_mean(data_test, 'fixed', W);
                    
                    ISC_av_pca(ik, irep) = mean(ISC(1:K, ik, irep));
                    ISC_av_test_pca(ik, irep) = mean(ISC_test(1:K, ik, irep));
                    
                end
                ik
            end
            
            % loading stuff precomputd by Stefan
            load('stefan/sims/mat/sim1a_results.mat','ISC_av','ISC_av_test');
            
            clf
            h0 = plot([min(Ks) max(Ks)], [1 1], '--', 'linewidth', 2, 'Markersize', 8, 'Color', 'k');
            hold on
            set(gca, 'ColorOrderIndex', 4)
            h1 = plot(Ks, [mean(ISC_av, 2)], 's-', 'linewidth', 2, 'Markersize', 8);
            set(h1, 'MarkerFaceColor', get(h1, 'Color'));
            patch([Ks'; fliplr(Ks)'], [mean(ISC_av, 2) + std(ISC_av, [], 2); flipud(mean(ISC_av, 2) - std(ISC_av, [], 2))], ...
                get(h1, 'Color'), 'edgecolor', 'none', 'facealpha', 0.1)
            h2 = plot(Ks, [mean(ISC_av_test, 2)], 's-', 'linewidth', 2, 'Markersize', 8);
            set(h2, 'MarkerFaceColor', get(h2, 'Color'));
            patch([Ks'; fliplr(Ks)'], [mean(ISC_av_test, 2) + std(ISC_av_test, [], 2); flipud(mean(ISC_av_test, 2) - std(ISC_av_test, [], 2))], ...
                get(h2, 'Color'), 'edgecolor', 'none', 'facealpha', 0.1)
            h3 = plot(Ks, [mean(ISC_av_pca, 2)], 's-', 'linewidth', 2, 'Markersize', 8);
            set(h3, 'MarkerFaceColor', get(h3, 'Color'));
            patch([Ks'; fliplr(Ks)'], [mean(ISC_av_pca, 2) + std(ISC_av_pca, [], 2); flipud(mean(ISC_av_pca, 2) - std(ISC_av_pca, [], 2))], ...
                get(h3, 'Color'), 'edgecolor', 'none', 'facealpha', 0.1)
            h4 = plot(Ks, [mean(ISC_av_test_pca, 2)], 's-', 'linewidth', 2, 'Markersize', 8);
            set(h4, 'MarkerFaceColor', get(h4, 'Color'));
            patch([Ks'; fliplr(Ks)'], [mean(ISC_av_test_pca, 2) + std(ISC_av_test_pca, [], 2); flipud(mean(ISC_av_test_pca, 2) - std(ISC_av_test_pca, [], 2))], ...
                get(h4, 'Color'), 'edgecolor', 'none', 'facealpha', 0.1)
            
            set(gca, 'fontsize', 20)
            grid on
            xlim([min(Ks) max(Ks)])
            ylim([0 1.05])
            xlabel('True K')
            ylabel('Mean ISC')
            title('IID')
            legend([h1 h2 h3 h4 h0], 'Train', 'Test',  'Train-PCA', 'Test-PCA', 'Ideal', 'Location', 'Southwest')
            ax = reshape(axis, 2, 2);
            text(ax(1, 1) - diff(ax(:, 1))*0.125, ax(1, 2) + diff(ax(:, 2))*1.025, 'A', 'fontsize', 30)
            export_fig('figures/simulations/sim1a_ISC_pca_mean', '-r300', '-a2','-transparent');
            
            % SNR
            db_snrs = [-40:10:40];
            snrs = 10.^(db_snrs/20)./(1+10.^(db_snrs/20)); %=0.5
            
            % K: number of correlated components
            K = 10;
            
            clear ISC ISC_av_pca ISC_av_test_pca
            for isnr = 1:length(snrs)
                snr = snrs(isnr);
                
                for irep = 1:Nrep
                    
                    [data, Y_true, A_true] = gen_data(2*T, D, N, ones(1, K), snr, 0, 1, 1, 0, 'Gauss');
                    
                    data_test = data(T+1:end, :, :);
                    data = data(1:T, :, :);
                    
                    [W, ISC(:,isnr,irep)] = pca_mean(data);
                    [~, ISC_test(:,isnr,irep)] = pca_mean(data_test, 'fixed', W);
                    
                    ISC_av_pca(isnr, irep) = mean(ISC(1:K, isnr, irep));
                    ISC_av_test_pca(isnr, irep) = mean(ISC_test(1:K, isnr, irep));
                    
                end
                isnr
            end
            
            % loading stuff precomputd by Stefan
            load('stefan/sims/mat/sim2a_results.mat','ISC_av','ISC_av_test');
            
            clf
            h0 = plot([min(db_snrs) max(db_snrs)], [1 1], '--', 'linewidth', 2, 'Markersize', 8, 'Color', 'k');
            hold on
            set(gca, 'ColorOrderIndex', 4)
            h1 = plot(db_snrs, [mean(ISC_av, 2)], 's-', 'linewidth', 2, 'Markersize', 8);
            set(h1, 'MarkerFaceColor', get(h1, 'Color'));
            patch([db_snrs'; fliplr(db_snrs)'], [mean(ISC_av, 2) + std(ISC_av, [], 2); flipud(mean(ISC_av, 2) - std(ISC_av, [], 2))], ...
                get(h1, 'Color'), 'edgecolor', 'none', 'facealpha', 0.1)
            h2 = plot(db_snrs, [mean(ISC_av_test, 2)], 's-', 'linewidth', 2, 'Markersize', 8);
            set(h2, 'MarkerFaceColor', get(h2, 'Color'));
            patch([db_snrs'; fliplr(db_snrs)'], [mean(ISC_av_test, 2) + std(ISC_av_test, [], 2); flipud(mean(ISC_av_test, 2) - std(ISC_av_test, [], 2))], ...
                get(h2, 'Color'), 'edgecolor', 'none', 'facealpha', 0.1)
            h3 = plot(db_snrs, [mean(ISC_av_pca, 2)], 's-', 'linewidth', 2, 'Markersize', 8);
            set(h3, 'MarkerFaceColor', get(h1, 'Color'));
            patch([db_snrs'; fliplr(db_snrs)'], [mean(ISC_av_pca, 2) + std(ISC_av_pca, [], 2); flipud(mean(ISC_av_pca, 2) - std(ISC_av_pca, [], 2))], ...
                get(h3, 'Color'), 'edgecolor', 'none', 'facealpha', 0.1)
            h4 = plot(db_snrs, [mean(ISC_av_test_pca, 2)], 's-', 'linewidth', 2, 'Markersize', 8);
            set(h4, 'MarkerFaceColor', get(h4, 'Color'));
            patch([db_snrs'; fliplr(db_snrs)'], [mean(ISC_av_test_pca, 2) + std(ISC_av_test_pca, [], 2); flipud(mean(ISC_av_test_pca, 2) - std(ISC_av_test_pca, [], 2))], ...
                get(h4, 'Color'), 'edgecolor', 'none', 'facealpha', 0.1)
            
            set(gca, 'fontsize', 20)
            grid on
            xlim([min(db_snrs) max(db_snrs)])
            ylim([0 1.05])
            xlabel('SNR [dB]')
            ylabel('Mean ISC')
            title('IID')
            legend([h1 h2 h3 h4 h0], 'Train', 'Test',  'Train-PCA', 'Test-PCA', 'Ideal', 'Location', 'Southeast')
            ax = reshape(axis, 2, 2);
            text(ax(1, 1) - diff(ax(:, 1))*0.125, ax(1, 2) + diff(ax(:, 2))*1.025, 'A', 'fontsize', 30)
            export_fig('figures/simulations/sim2a_ISC_pca_mean', '-r300', '-a2','-transparent');
            
        case 'NBDT Reviewer 1 hypothetical example'
            
            set1=randn(20,1);
            set2=randn(20,1);
            q1 = [set1, set2]; % questions 1, parent and child
            q2 = [-set2, -set1]; % questions 2, parent and child
            M =  [-0.7071    0.7071; 0.7071    0.7071];
            
            X(:,:,1) = [q1(:,1) q2(:,1)]; % parents
            X(:,:,2) = [q1(:,2) q2(:,2)]; % children
            
            D=2;
            
            ISCx=corrcoef([X(:,:,1),X(:,:,2)]); ISCx=diag(ISCx(1:D,D+1:2*D));;
            
            % intersubject correlation optimized by CorrCA
            %             [W,ISCy,Y]=corrca(X,'version',2);
            [K, p, W, ISCy, Y, A] = est_K_shift(X, 1000, 0);
            
            % show all time courses
            clf, clear h;
            for i=1:D, h(i)    =subplot(2*D,2,i*2-1); plot(X(:,i,1),'v-'); hold on; plot(X(:,i,2),'^-'); text(0.6,0.9,['ISC=' num2str(ISCx(i),1)],'units','normalized','backgroundcolor',[1 1 1]); ylabel(['X_' num2str(i)]); end; xlabel('time (samples)')
            for i=1:D, h(i+D+1)=subplot(2*D,2,i*2  ); plot(Y(:,i,1),'v-'); hold on; plot(Y(:,i,2),'^-'); text(0.6,0.9,['ISC=' num2str(ISCy(i),1) ', p=' num2str(p(i),1)],'units','normalized','backgroundcolor',[1 1 1]); ylabel(['Y_' num2str(i)]); end; xlabel('time (samples)')
            subplot(2*D,2,1); title('Original data, X');
            subplot(2*D,2,2); title('Correlated Components, Y=W^T X');
            hlegend=legend('Subject 1','Subject 2'); hlegend.Position([1 2])=[0.45 0.476];
            
            % show first two dimensions that have been rotated as scatter plot
            h(D+1)=subplot(2,2,3);
            
            plot(X(:,1,1),X(:,2,1),'v',X(:,1,2),X(:,2,2),'^'); hold on
            plot([X(:,1,1) X(:,1,2)]',[X(:,2,1) X(:,2,2)]','color',[0.8 0.8 0.8]); hold off
            axis equal;
            xlabel('X_1'); ylabel('X_2')
            arrow([0 0], M(1,:), 'color','r');
            arrow([0 0], M(2,:), 'color','b');
            
            
            h(2*D+2)=subplot(2,2,4);
            plot(Y(:,1,1),Y(:,2,1),'v',Y(:,1,2),Y(:,2,2),'^'); hold on
            plot([Y(:,1,1) Y(:,1,2)]',[Y(:,2,1) Y(:,2,2)]','color',[0.8 0.8 0.8]);
            axis equal; hold off
            xlabel('Y_1'); ylabel('Y_2')
            
            sublabel(h);
            
        case 'NBDT Reviewer 2 question'
            
            % This is Lucas' simulation of what happens in the common
            % source is delayed. Stefan did not agree that this is a good
            % example and saw more merit in the topic. In the final
            % manuscript/response Stefan pravailed in the answer to this
            % question from the rievwer. So this is here just so we dont
            % forget the conversation. 
            
            T=100;
            t=(1:T)'/T;
            
            alpha1=pi/4;alpha2=pi/3;
            M = [sin(alpha1) cos(alpha1);cos(alpha2) -sin(alpha2)]*diag([1 1]);
            
            a=0.5;b=0.5;c=-0.5; d=0.5;
            for l=1:30
                s = abs(sin(2*pi*2*t)*(rand(1)-0.5)+ cos(2*pi*2*t)*(rand(1)-0.5));
                n = randn(T,1);
                X(:,:,l) = [a*s + c*n, b*s + d*n];
            end
            [W,ISC,Y,A] = corrca(X);
            
            plot(mean(Y,3))
            
    end
    
    str = test; str(str==' ')='-';
    if savefigure, saveas(gcf,[figurepath str],'png'); saveas(gcf,[figurepath str],'epsc'); end
    
end

