#include "HistoHandler/HistoHandler.h"

HistoHandler& HistoHandler::getInstance( )
{
    static HistoHandler hists;
    return hists;
}
