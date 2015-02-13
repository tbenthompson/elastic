#include "controller.h"

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage is 'solve filename'" << std::endl;
        return 1;
    }
    
    auto filename = argv[1];
    solve<3>(filename);
}
