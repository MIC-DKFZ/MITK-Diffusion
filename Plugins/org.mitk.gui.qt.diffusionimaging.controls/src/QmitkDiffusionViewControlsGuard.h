/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef QmitkDiffusionViewControlsGuard_h
#define QmitkDiffusionViewControlsGuard_h

#include <QObject>
#include <QWidget>

#include <functional>

/**
 * \brief Disconnects a view from its widgets before they are deleted.
 *
 * BlueBerry deletes the widgets of a closed view before the view itself. MITK's
 * QmitkDataStorageComboBox empties itself in its destructor and so emits currentIndexChanged and
 * OnSelectionChanged; a view slot connected to them then runs while the view's other widgets
 * are being deleted, and touching one that is already gone crashes the Workbench (closing the
 * view, or quitting, with data loaded). Which widgets are gone depends on their order in the
 * .ui file.
 */
namespace QmitkDiffusionViewControlsGuard
{
  /**
   * \brief Disconnects every widget below parent from view when parent deletes its children.
   *
   * Call it in CreateQtPartControl() with its parent before setupUi(parent): the guard is then
   * the first child of parent, and Qt deletes children in the order they were created, so the
   * guard goes while all widgets still exist. A view whose timers or worker threads touch its
   * widgets passes beforeDestruction, which runs then too and stops them.
   */
  inline void Install(QWidget *parent, QObject *view, std::function<void()> beforeDestruction = {})
  {
    auto *guard = new QObject(parent);
    QObject::connect(guard, &QObject::destroyed, view, [parent, view, beforeDestruction]() {
      if (beforeDestruction)
        beforeDestruction();
      // While parent deletes its children, the entry of the child being deleted (the guard) is
      // null, so the list is walked by hand.
      for (auto *child : parent->children())
      {
        if (child == nullptr)
          continue;
        child->disconnect(view);
        for (auto *descendant : child->findChildren<QObject *>())
          descendant->disconnect(view);
      }
    });
  }
}

#endif
