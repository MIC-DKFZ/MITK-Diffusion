/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef QmitkReconstructionResults_h
#define QmitkReconstructionResults_h

#include <mitkDataNode.h>
#include <mitkLevelWindowProperty.h>
#include <mitkLog.h>
#include <mitkLookupTableProperty.h>

#include <itkMacro.h>

#include <QMessageBox>

#include <exception>
#include <string>

/** What the views of this plugin share: result nodes and the report of a failed action. */
namespace QmitkReconstructionResults
{
  /** A node for data computed from parent, named after it with the suffix; the caller adds it. */
  inline mitk::DataNode::Pointer ResultNode(mitk::BaseData *data, const mitk::DataNode *parent, const std::string &suffix)
  {
    auto node = mitk::DataNode::New();
    node->SetData(data);
    node->SetName(parent->GetName() + suffix);
    return node;
  }

  /** Jet colours, transparent at zero, and the given level window, as for scalar index maps. */
  inline void ShowAsIndexMap(mitk::DataNode *node, double level, double window)
  {
    auto lut = mitk::LookupTable::New();
    lut->SetType(mitk::LookupTable::JET_TRANSPARENT);
    auto lutProperty = mitk::LookupTableProperty::New();
    lutProperty->SetLookupTable(lut);
    node->SetProperty("LookupTable", lutProperty);

    mitk::LevelWindow levelWindow;
    levelWindow.SetLevelWindow(level, window);
    node->SetProperty("levelwindow", mitk::LevelWindowProperty::New(levelWindow));
  }

  /**
   * Runs the action of a button. An exception is shown in a message box with the title: through Qt
   * it would end in MITK's dialog that offers to quit the application.
   */
  template <typename F>
  void Run(const char *title, F &&action)
  {
    try
    {
      action();
    }
    catch (const itk::ExceptionObject &e)
    {
      MITK_ERROR << e;
      QMessageBox::information(nullptr, title, e.GetDescription());
    }
    catch (const std::exception &e)
    {
      MITK_ERROR << e.what();
      QMessageBox::information(nullptr, title, e.what());
    }
  }
}

#endif
