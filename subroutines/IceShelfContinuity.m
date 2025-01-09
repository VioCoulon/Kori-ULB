function [H, Hn] = IceShelfContinuity(ctr, row, col, H, Hn, glMASK)

    for r=1:length(row)

        i = row(r);
        j = col(r);
        
        im = i;
        jm = j;
        ip = i;
        jp = j;

        for k=1:ctr.imax

            im = im - k;
            ip = ip + k;
            jp = jp + k;
            jm = jm - k;

            % Ensure indices within limits.
            if (im < 2) || (jm < 2) || (ip > ctr.imax-1) || (jp > ctr.jmax-1)
                break;
            end

            if glMASK(im,j)==4 || glMASK(im,j)==3 %&& glMASK_old(im,j)==4 instead??

                H(i,j)=H(im,j); % Necessary to add this?
                Hn(i,j)=Hn(im,j);
                %Hn(i,j)=2*Hn(im,j) - Hn(im-1,j);
                %H(i,j)=Hn(i,j);

                break; 
            end

            if glMASK(ip,j)==4 || glMASK(ip,j)==3 %&& glMASK(ip,j)==4

                H(i,j)=H(ip,j);
                Hn(i,j)=Hn(ip,j);
                %Hn(i,j)=2*Hn(ip,j) - Hn(ip+1,j);
                %H(i,j)=Hn(i,j);

                break;
            end

            if glMASK(i,jm)==4 || glMASK(i,jm)==3 %&& glMASK(i,jm)==4

                H(i,j)=H(i,jm);
                Hn(i,j)=Hn(i,jm);
                %Hn(i,j)=2*Hn(i,jm) - Hn(i,jm-1);
                %H(i,j)=Hn(i,j);

                break;
            end
            
            if glMASK(i,jp)==4 || glMASK(i,jp)==3 %&& glMASK(i,jp)==4

                H(i,j)=H(i,jp);
                Hn(i,j)=Hn(i,jp);
                %Hn(i,j)=2*Hn(i,jp) - Hn(i,jp+1);
                %H(i,j)=Hn(i,j);

                break;
            end

            if glMASK(im,jm)==4 || glMASK(im,jm)==3 %&& glMASK(im,jm)==4                    

                H(i,j)=H(im,jm);
                Hn(i,j)=Hn(im,jm);
                %Hn(i,j)=2*Hn(im,jm) - Hn(im-1,jm-1);
                %H(i,j)=Hn(i,j);
                
                break;
            end

            if glMASK(im,jp)==4 || glMASK(im,jp)==3 %&& glMASK(im,jp)==4

                H(i,j)=H(im,jp);
                Hn(i,j)=Hn(im,jp);
                %Hn(i,j)=2*Hn(im,jp) - Hn(im-1,jp+1);
                %H(i,j)=Hn(i,j);
                %fprintf('\n im, jp')

                break;
            end

            if glMASK(ip,jm)==4 || glMASK(ip,jm)==3 %&& glMASK(ip,jm)==4

                H(i,j)=H(ip,jm);
                Hn(i,j)=Hn(ip,jm);
                %Hn(i,j)=2*Hn(ip,jm) - Hn(ip+1,jm-1);
                %H(i,j)=Hn(i,j);
                %fprintf('\n ip, jm')

                break;
            end

            if glMASK(ip,jp)==4 || glMASK(ip,jp)==3 %&& glMASK(ip,jp)==4

                H(i,j)=H(ip,jp);
                Hn(i,j)=Hn(ip,jp);
                %Hn(i,j)=2*Hn(ip,jp) - Hn(ip+1,jp+1);
                %H(i,j)=Hn(i,j);
                %fprintf('\n ip,ip')

                break;
            end
        end

    end

end
