// Links against the Core layer alone: IntStr's append is defined in the library.
#include <igor/Core.h>

int main()
{
    igor::core::IntStr sequence;
    sequence.append(3);
    return sequence.size() == 1 ? 0 : 1;
}
