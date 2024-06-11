#include <Exception.h>
#include <MathematicalOperations.h>
#include <MathematicalOperationsSTL.h>


MathematicalOperations* MathematicalOperations::getInstance(MathImplementaionType type)
{
    if(MathImplementaionType::BASE == type)
    {
        return MathematicalOperations::getInstanceBase();
    }
    else if(MathImplementaionType::STD == type)
    {
        return MathematicalOperationsSTL::getInstance();
    }

    throw RuntimeException("Unknown MathematicalOperations type!!");
}