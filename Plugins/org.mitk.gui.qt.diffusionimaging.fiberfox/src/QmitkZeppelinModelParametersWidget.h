/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
#ifndef _QMITKZeppelinModelParametersWidget_H_INCLUDED
#define _QMITKZeppelinModelParametersWidget_H_INCLUDED

//QT headers
#include <QWidget>
#include <QString>
#include "ui_QmitkZeppelinModelParametersWidgetControls.h"
#include <org_mitk_gui_qt_diffusionimaging_fiberfox_Export.h>

/** @brief
  */
class DIFFUSIONIMAGING_FIBERFOX_EXPORT QmitkZeppelinModelParametersWidget : public QWidget
{
    //this is needed for all Qt objects that should have a MOC object (everything that derives from QObject)
    Q_OBJECT

public:

    static const std::string VIEW_ID;

    QmitkZeppelinModelParametersWidget (QWidget* parent = nullptr, Qt::WindowFlags f = Qt::WindowFlags());
    virtual ~QmitkZeppelinModelParametersWidget();

    virtual void CreateQtPartControl(QWidget *parent);

    void SetT1(double t1){ m_Controls->m_T1box->setValue(t1); }
    void SetT2(double t2){ m_Controls->m_T2box->setValue(t2); }
    void SetD1(double d1){ m_Controls->m_D1box->setValue(d1); }
    void SetD2(double d2){ m_Controls->m_D2box->setValue(d2); }

    double GetD1(){ return m_Controls->m_D1box->value(); }
    double GetD2(){ return m_Controls->m_D2box->value(); }
    unsigned int GetT2(){ return m_Controls->m_T2box->value(); }
    unsigned int GetT1(){ return m_Controls->m_T1box->value(); }

public slots:

    void DChanged( double value );

protected:
    // member variables
    Ui::QmitkZeppelinModelParametersWidgetControls* m_Controls;
    void UpdateUi();

private:

};

#endif // _QMITKZeppelinModelParametersWidget_H_INCLUDED

