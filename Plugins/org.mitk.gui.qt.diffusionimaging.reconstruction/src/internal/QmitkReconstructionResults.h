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
#include <mitkDataStorage.h>
#include <mitkLevelWindowProperty.h>
#include <mitkLog.h>
#include <mitkLookupTableProperty.h>

#include <itkMacro.h>

#include <QFutureWatcher>
#include <QMessageBox>
#include <QtConcurrentRun>
#include <QWidget>

#include <exception>
#include <string>
#include <type_traits>
#include <utility>

/** What the views of this plugin share: result nodes, the report of a failed action and running it
 * in the background. */
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

  /** Adds the node below parent, or at the top level if parent was removed meanwhile. */
  inline void AddResult(mitk::DataStorage *dataStorage, mitk::DataNode *node, mitk::DataNode *parent)
  {
    dataStorage->Add(node, parent != nullptr && dataStorage->Exists(parent) ? parent : nullptr);
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

  /**
   * Runs work in a worker thread, so that the GUI stays responsive, and then done(result) in the GUI
   * thread; an exception of either is reported as by Run. A view runs one action at a time: busy is
   * true from now until done has run, and update() is called whenever busy changes, e.g. to disable
   * the start buttons. work must neither touch widgets nor the data storage, and holds its own
   * (const) smart pointers to its input, since nodes may be removed while it runs; reading through
   * a const image takes a read accessor, which does not block rendering. The watcher is a child of
   * context, a widget of the view: BlueBerry deletes the widgets of a closed view before the view,
   * so once they are gone, nothing else is called and the result is dropped. Nothing
   * cancels work (the filters cannot be interrupted): it keeps its pool thread busy until it ends,
   * also after the view was closed, and quitting the application waits for it, because Qt waits
   * for its global thread pool on exit.
   */
  template <typename Update, typename Work, typename Done>
  void RunInBackground(QWidget *context, bool &busy, Update update, const char *title, Work work, Done done)
  {
    using Result = std::invoke_result_t<Work &>;
    struct Outcome
    {
      Result result{};
      std::exception_ptr error;
    };

    busy = true;
    update();
    auto watcher = new QFutureWatcher<Outcome>(context);
    QObject::connect(watcher, &QFutureWatcherBase::finished, context,
                     [watcher, &busy, update, title, done = std::move(done)]() mutable {
                       // Taken, not copied: the future would keep the result alive.
                       Outcome outcome = watcher->future().takeResult();
                       watcher->deleteLater();
                       Run(title, [&]() {
                         if (outcome.error)
                           std::rethrow_exception(outcome.error);
                         done(std::move(outcome.result));
                       });
                       busy = false;
                       update();
                     });
    watcher->setFuture(QtConcurrent::run([work = std::move(work)]() mutable {
      Outcome outcome;
      try
      {
        outcome.result = work();
      }
      catch (...)
      {
        outcome.error = std::current_exception();
      }
      return outcome;
    }));
  }
}

#endif
