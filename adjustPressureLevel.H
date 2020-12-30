// Compute the full pressure and adjust the pressure level.
Boussinesq.updateBackgroundPressure();
p = p_rgh + Boussinesq.backgroundPressure();

if (p_rghNeedRef)
{
    Info<< "p_rgh[pRefCell] (before adjustment): " 
    	<< getRefCellValue(p_rgh, pRefCell) << endl;
    
    dimensionedScalar pAdjust(
    	"pAdjust",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell));

    p     += pAdjust;
    p_rgh += pAdjust;

    Info<< "p_rgh[pRefCell] (after  adjustment): " 
    	<< getRefCellValue(p_rgh, pRefCell) << endl;
}
