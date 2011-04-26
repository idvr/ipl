/*
 * IPL (Image Processing Library)
 *
 * Copyright (C) 2011 Andrea Gagliardi La Gala.
 * Email: andrea@ameliemedia.com
 * Web: http://www.ameliemedia.com/
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * IPL depends on and extends:
 *   - Leptonica library by Dan S. Bloomberg (http://www.leptonica.com/)
 *
 * and (optionally) integrates:
 *   - FFTW library by Matteo Frigo and Steven G. Johnson (http://www.fftw.org/)
 *   - fuzzylite C++ library by Juan Rada-Vilela (http://code.google.com/p/fuzzy-lite/)
 *
 * to achieve fast and efficient image filtering and manipulation within
 * the frequency domain (using the Discrete Fourier Transform and its inverse)
 * and by using Fuzzy Sets and Logic.
 *
 */

/*
 *  pixafunc.c
 *
 *   Filters
 *       PIX      *pixInnerSelectBySize()
 */

#include "ipl.h"

/*!
 *  pixInnerSelectBySize()
 *
 *      Input:  pixs (1 bpp)
 *              width, height (threshold dimensions)
 *              connectivity (4 or 8)
 *              type (L_SELECT_WIDTH, L_SELECT_HEIGHT,
 *                    L_SELECT_IF_EITHER, L_SELECT_IF_BOTH)
 *              relation (L_SELECT_IF_LT, L_SELECT_IF_GT,
 *                        L_SELECT_IF_LTE, L_SELECT_IF_GTE)
 *              &changed (<optional return> 1 if changed; 0 otherwise)
 *      Return: filtered pixd, or null on error
 *
 *  Notes:
 *      (1) The args specify constraints on the size of the
 *          components that are kept. Components of any size
 *          that touch the pix borders are always removed.
 *      (2) If unchanged, returns a copy of pixs.  Otherwise,
 *          returns a new pix with the filtered components.
 *      (3) If the selection type is L_SELECT_WIDTH, the input
 *          height is ignored, and v.v.
 *      (4) To keep small components, use relation = L_SELECT_IF_LT or
 *          L_SELECT_IF_LTE.
 *          To keep large components, use relation = L_SELECT_IF_GT or
 *          L_SELECT_IF_GTE.
 */
PIX *
pixInnerSelectBySize(PIX      *pixs,
					 l_int32   width,
					 l_int32   height,
					 l_int32   connectivity,
					 l_int32   type,
					 l_int32   relation,
					 l_int32  *pchanged)
{
	l_int32  w, h, x, y, i, n, empty, changed, count;
	BOXA    *boxa;
	PIX     *pixd;
	PIXA    *pixas, *pixad;
	
    PROCNAME("pixInnerSelectBySize");
	
    if (!pixs)
        return (PIX *)ERROR_PTR("pixs not defined", procName, NULL);
    if (connectivity != 4 && connectivity != 8)
        return (PIX *)ERROR_PTR("connectivity not 4 or 8", procName, NULL);
    if (type != L_SELECT_WIDTH && type != L_SELECT_HEIGHT &&
        type != L_SELECT_IF_EITHER && type != L_SELECT_IF_BOTH)
        return (PIX *)ERROR_PTR("invalid type", procName, NULL);
    if (relation != L_SELECT_IF_LT && relation != L_SELECT_IF_GT &&
        relation != L_SELECT_IF_LTE && relation != L_SELECT_IF_GTE)
        return (PIX *)ERROR_PTR("invalid relation", procName, NULL);
    if (pchanged) *pchanged = FALSE;
    
	/* Check if any components exist */
    pixZero(pixs, &empty);
    if (empty)
        return pixCopy(NULL, pixs);
	
	/* Identify and select the components */
    boxa = pixConnComp(pixs, &pixas, connectivity); 
    pixad = pixaSelectBySize(pixas, width, height, type, relation, &changed);
	
	/* Remove all components attached to the border */
	n = pixaGetCount(pixad);
	for (i = 0; i < n; i++) {
		pixaGetBoxGeometry(pixad, i, &x, &y, &w, &h);
		if (x <= 0 || y <= 0 || x + w >= pixs->w || y + h >= pixs->h) {
			pixaRemovePix(pixad, i);
			n--;
		}
	}
    boxaDestroy(&boxa);
    pixaDestroy(&pixas);
	
	/* Render the result */
    if (!changed) {
        pixaDestroy(&pixad);
        return pixCopy(NULL, pixs);
    }
    else {
        if (pchanged) *pchanged = TRUE;
        pixGetDimensions(pixs, &w, &h, NULL);
        count = pixaGetCount(pixad);
        if (count == 0)  /* return empty pix */
            pixd = pixCreateTemplate(pixs);
        else {
            pixd = pixaDisplay(pixad, w, h);
            pixCopyResolution(pixd, pixs);
            pixCopyColormap(pixd, pixs);
            pixCopyText(pixd, pixs);
            pixCopyInputFormat(pixd, pixs);
        }
        pixaDestroy(&pixad);
        return pixd;
    }
}