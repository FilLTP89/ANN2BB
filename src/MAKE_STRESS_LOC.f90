subroutine MAKE_STRESS_LOC(lambda,mu,dE,dS)
    implicit none
            dS(1)=(lambda+2.0d0*mu)*dE(1)+lambda*dE(2)
            dS(2)=(lambda+2.0d0*mu)*dE(2)+lambda*dE(1)  
            dS(3)=mu*dE(3)
    return
end subroutine MAKE_STRESS_LOC
