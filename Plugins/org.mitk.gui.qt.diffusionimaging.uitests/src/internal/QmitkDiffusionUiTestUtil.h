/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef QmitkDiffusionUiTestUtil_h
#define QmitkDiffusionUiTestUtil_h

#include <berryIWorkbenchPage.h>

#include <mitkDataStorage.h>
#include <mitkLogBackendBase.h>

#include <cppunit/TestAssert.h>

#include <QColor>
#include <QList>
#include <QMessageBox>
#include <QStringList>
#include <QWidget>

#include <deque>
#include <functional>
#include <mutex>
#include <string>
#include <vector>

class QTimer;

/**
 * \brief Helpers of the diffusion UI tests, which run in the GUI thread of an open Workbench.
 *
 * The tests use the views as a user does: they open a view, find its widgets by the object names
 * of its .ui file, select data through the workbench selection and click. Nothing here links a
 * view's plug-in.
 */
namespace QmitkDiffusionUiTest
{
  berry::IWorkbenchPage::Pointer GetPage();

  /**
   * \brief Opens the view and returns the root widget of its controls.
   * \param controlsName Object name of the root widget of the view's .ui form; setupUi() gives
   * the view's parent widget this name.
   */
  QWidget *ShowView(const QString &viewId, const QString &controlsName);

  /** Closes the view if it is open; a closed view's controls are deleted. */
  void HideView(const QString &viewId);

  bool IsViewOpen(const QString &viewId);

  template <class T>
  T *GetChild(QWidget *root, const char *name)
  {
    auto *child = root->findChild<T *>(name);
    CPPUNIT_ASSERT_MESSAGE(std::string("No widget ") + name, child != nullptr);
    return child;
  }

  mitk::DataStorage::Pointer GetDataStorage();

  /** Removes every node, as closing a project does. */
  void ClearDataStorage();

  /** Loads a file of the diffusion test data into the data storage and returns its node. */
  mitk::DataNode::Pointer Load(const std::string &testDataFile);

  mitk::DataNode::Pointer AddNode(mitk::BaseData *data, const std::string &name, mitk::DataNode *parent = nullptr);

  /** The nodes whose data has this class name (GetNameOfClass()), in data storage order. */
  std::vector<mitk::DataNode::Pointer> GetNodesOfClass(const std::string &className);

  /**
   * \brief Chooses the node (nullptr: none) in a QmitkDataStorageComboBox of the view and checks
   * that the box holds it (it lists only nodes its predicate accepts).
   */
  void ChooseNode(QWidget *root, const char *boxName, mitk::DataNode *node);

  /** Selects the nodes as the Data Manager does (workbench selection service). */
  void Select(const QList<mitk::DataNode::Pointer> &nodes);

  /** Fits the render windows to the data, as the Workbench does after loading a file. */
  void InitializeViews();

  /** Moves the crosshair (the selected position of the render windows) to the world point. */
  void SetCrosshair(const mitk::Point3D &position);

  /** Processes events until the condition holds; false after the timeout. */
  bool WaitFor(const std::function<bool()> &condition, int timeoutMs = 120000);

  /** Processes events for the given time. */
  void Wait(int ms);

  /** A new folder below the system temporary folder. */
  std::string CreateTempDir();

  /**
   * \brief Answers the modal dialogs the views open.
   *
   * A test cannot answer a modal dialog itself: the click that opens it returns only when the
   * dialog is closed. A timer in the dialog's event loop closes it with the next queued answer, or
   * cancels it (message boxes: their escape button). The texts of all message boxes are kept.
   */
  class DialogHandler
  {
  public:
    static DialogHandler &Instance();

    void Start();
    void Stop();

    void AnswerMessageBox(QMessageBox::StandardButton button);
    void AnswerFileDialog(const QString &path);
    void AnswerColorDialog(const QColor &color);

    /** Texts (title and text) of the message boxes closed since the last call. */
    QStringList TakeMessages();

    /** The same, without clearing them. */
    QStringList PeekMessages() const;

  private:
    DialogHandler() = default;
    void HandleModalWidget();

    QTimer *m_Timer = nullptr;
    std::deque<std::function<bool(QWidget *)>> m_Answers;
    QStringList m_Messages;
    bool m_Busy = false;
  };

  /** \brief Collects the MITK log messages while it exists. */
  class LogCapture : public mitk::LogBackendBase
  {
  public:
    LogCapture();
    ~LogCapture() override;

    void ProcessMessage(const mitk::LogMessage &message) override;
    OutputType GetOutputType() const override;

    bool Contains(const std::string &text) const;
    std::vector<std::string> GetMessages() const;

    /** True if an error or fatal message was logged. */
    bool ContainsError() const;

  private:
    mutable std::mutex m_Mutex;
    std::vector<std::string> m_Messages;
    std::vector<std::string> m_Errors;
  };
}

#endif
