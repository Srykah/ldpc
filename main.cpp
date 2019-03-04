#include <iostream>
#include "fonctions.hpp"

int main()
{
    Matrice test({
        { 0, 5, 7, 9 },
        { 0, 3, 4, 10 },
        { 1, 4, 6, 8 },
        { 2, 5, 10, 11 },
        { 2, 6, 7, 11 },
        { 1, 4, 8, 10 },
        { 0, 3, 6, 9 },
        { 1, 5, 7, 9 },
        { 2, 3, 8, 11 }
    });

    Matrice test2({
        { 1, 2, 4, 5 },
        { 0, 1, 3, 4 },
        { 0, 2, 3, 5 }
    });

    // gauss_jordan(test);
    // std::cout << "\ntest :\n" << test << "\ntest2 :\n" << test2 << std::endl;

    Matrice test3 = gallager(12, 4, 9);

    std::cout << "\ntest3 : \n" << test3 << std::endl;

    return 0;
}
