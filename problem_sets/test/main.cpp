#include <iostream>

#include "fact.h"

using std::cout;
using std::endl;

int main() {

    cout << "Hello World!" << endl;

    for (int i = 1; i < 10; i++) {

        cout << i << " : " << fact(i) << endl;

    }

    return 0;

}
