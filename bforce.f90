SUBROUTINE BFORCE
!### subroutine to calculate body force per unit volume
!###
   USE DATA
   IMPLICIT NONE
   
   integer::It,Iteration,Num_Point,Num_Elem
   integer::J_start,J_end
   integer::coord_1,coord_2,coord_3,coord_4
   real(dp),allocatable::X_s(:),Y_s(:),Fx_s(:),Fy_s(:)
   integer,allocatable::conn(:,:)
   real(dp)::x_min,x_max,y_min,y_max
   real(dp)::Fx_bot,Fx_top,Fy_bot,Fy_top
   real(dp)::force_max
   character(len=30)::filename

   if(NT>2) goto 99

   do Iteration=1,2
  
      if(Iteration==1) then
         filename='suzen_force_up.dat'
         J_start=JB;J_end=NJ
      else
         filename='suzen_force_bot.dat'
         J_start=1;J_end=JB-1
      end if

      open(unit=50,file=filename)
      
      read(50,*) Num_Point,Num_Elem
      allocate(X_s(Num_Point),Y_s(Num_Point),Fx_s(Num_Point),Fy_s(Num_Point))
      allocate(conn(Num_Elem,4))
      
      do It=1,Num_Point
            read(50,*) X_s(It),Y_s(It),Fx_s(It),Fy_s(It)
      end do
      
      do It=1,Num_Elem
         read(50,*) conn(It,1),conn(It,2),conn(It,3),conn(It,4)
      end do
      
      close(50)


      ! determine the maximum and minimum of the flow mesh
      x_min=10000
      x_max=0
      y_min=10000
      y_max=0
      do It=1,Num_Point
         if(x_s(It)<x_min) then 
             x_min=x_s(It)
         end if
         if(x_s(It)>x_max) then
             x_max=x_s(It)
         end if
         if(y_s(It)<y_min) then
             y_min=y_s(It)
         end if
         if(y_s(It)>y_max) then
             y_max=y_s(It)
         end if
      end do
    
   
      do I=1,NI
           if(xu(I)>=x_min.and.xu(I)<=x_max) then
     
             do J=J_start,J_end
               if(y(J)>=y_min.and.y(J)<=y_max) then
     
                 do It=1,Num_Elem

                    if(Iteration==1) then
                       coord_1=conn(It,1);coord_2=conn(It,2);coord_3=conn(It,3);coord_4=conn(It,4)
                    else
                       coord_1=conn(It,4);coord_2=conn(It,3);coord_3=conn(It,2);coord_4=conn(It,1)
                    end if

                    if(xu(I)>=x_s(coord_1).and.xu(I)<=x_s(coord_2).and.y(J)>=y_s(coord_4).and.y(J)<=y_s(coord_1)) then
                      
                       Fx_top=Fx_s(coord_1)+(xu(I)-x_s(coord_1))/(x_s(coord_2)-x_s(coord_1))&
                              &*(Fx_s(coord_2)-Fx_s(coord_1))
                       Fx_bot=Fx_s(coord_4)+(xu(I)-x_s(coord_4))/(x_s(coord_3)-x_s(coord_4))&
                              &*(Fx_s(coord_3)-Fx_s(coord_4))
                       Force_x(I,J)=Fx_bot+(y(J)-y_s(coord_4))/(y_s(coord_1)-y_s(coord_4))*(Fx_top-Fx_bot)
     
                    end if
                  end do
     
               end if
             end do
     
           end if
    
 
           if(x(I)>=x_min.and.x(I)<=x_max) then
     
             do J=J_start,J_end
               if(yv(J)>=y_min.and.yv(J)<=y_max) then
     
                 do It=1,Num_Elem

                    if(Iteration==1) then
                       coord_1=conn(It,1);coord_2=conn(It,2);coord_3=conn(It,3);coord_4=conn(It,4)
                    else
                       coord_1=conn(It,4);coord_2=conn(It,3);coord_3=conn(It,2);coord_4=conn(It,1)
                    end if

                    if(x(I)>=x_s(coord_1).and.x(I)<=x_s(coord_2).and.yv(J)>=y_s(coord_4).and.yv(J)<=y_s(coord_1)) then
                      
                       Fy_top=Fy_s(coord_1)+(x(I)-x_s(coord_1))/(x_s(coord_2)-x_s(coord_1))&
                              &*(Fy_s(coord_2)-Fy_s(coord_1))
                       Fy_bot=Fy_s(coord_4)+(x(I)-x_s(coord_4))/(x_s(coord_3)-x_s(coord_4))&
                              &*(Fy_s(coord_3)-Fy_s(coord_4))
                       Force_y(I,J)=Fy_bot+(yv(J)-y_s(coord_4))/(y_s(coord_1)-y_s(coord_4))*(Fy_top-Fy_bot)
     
                    end if
                  end do
     
               end if
             end do
     
           end if
     
       end do

       deallocate(x_s,y_s,Fx_s,Fy_s)
       deallocate(conn)   
   
   end do


   
   
99     do I=1,NI
          do J=1,NJ
             FBX(I,J)=FBSCALE*Force_x(I,J)*(1.0D0+sign(1.0D0,(sin(2*pi*FBFHZ*time+(0.5D0-tau)*pi)-sin((0.5D0-tau)*pi))))/2.0D0
             FBY(I,J)=FBSCALE*Force_y(I,J)*(1.0D0+sign(1.0D0,(sin(2*pi*FBFHZ*time+(0.5D0-tau)*pi)-sin((0.5D0-tau)*pi))))/2.0D0
!             FBX(I,J)=FBSCALE*Force_x(I,J)*(1.0D0+sign(1.0D0,sin(2*pi*FBFHZ*time)))/2.0D0
!             FBY(I,J)=FBSCALE*Force_y(I,J)*(1.0D0+sign(1.0D0,sin(2*pi*FBFHZ*time)))/2.0D0
          end do
       end do

     
 

END SUBROUTINE BFORCE






















