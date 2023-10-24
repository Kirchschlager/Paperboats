module functs
  use omp_lib
  use datatype
  use variable
  
  implicit none  
  
    !--------------------------------------------------------------------------!
    PUBLIC :: fun, cross_product, f_norm
    !--------------------------------------------------------------------------!
CONTAINS   

FUNCTION fun (x) result(fun_res) ! #s1.3.1.1
    use omp_lib
    use datatype
    use variable

    implicit none
    real(kind=r2), intent(in):: x
    real(kind=r2)            :: fun_res   
    
    fun_res = exp(x) - (m_electron/(mu_chem * m_amu))**0.5 * (1.0 - x)
 
END FUNCTION fun  ! #s1.3.1.1


function cross_product(x, y)  
    real(kind=r2), dimension(3), intent(in) :: x, y
    real(kind=r2), dimension(3)             :: cross_product     
 
    cross_product(1) = x(2)*y(3) - x(3)*y(2)
    cross_product(2) = x(3)*y(1) - x(1)*y(3)
    cross_product(3) = x(1)*y(2) - x(2)*y(1)
    
end function cross_product

function f_norm(x)  
    real(kind=r2), dimension(3), intent(in) :: x 
    real(kind=r2)                           :: f_norm    
 
    f_norm = sqrt(sum(x**2.0))
    
end function f_norm


end module functs


