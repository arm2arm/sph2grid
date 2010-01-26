/* 
 * File:   Render.h
 * Author: arm2arm
 *
 * Created on January 26, 2010, 5:44 PM
 */

#ifndef _RENDER_H
#define	_RENDER_H

class CRender {
public:
    CRender();
    CRender(const CRender& orig);
    virtual ~CRender();
    void DoRender(float ***vol3d, int GRID=128);
private:

};

#endif	/* _RENDER_H */

