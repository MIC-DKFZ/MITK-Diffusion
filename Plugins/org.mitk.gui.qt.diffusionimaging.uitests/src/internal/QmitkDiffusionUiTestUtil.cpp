/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "QmitkDiffusionUiTestUtil.h"

#include "QmitkDiffusionUiTestSelectionView.h"

#include <berryIWorkbenchWindow.h>
#include <berryPlatformUI.h>

#include <mitkCoreServices.h>
#include <mitkDiffusionTestData.h>
#include <mitkIDataStorageService.h>
#include <mitkIOUtil.h>
#include <mitkLog.h>
#include <mitkRenderingManager.h>
#include <mitkWorkbenchUtil.h>

#include <QmitkDataStorageComboBoxWithSelectNone.h>

#include <QAbstractButton>
#include <QApplication>
#include <QColorDialog>
#include <QFileDialog>
#include <QFileInfo>
#include <QLineEdit>
#include <QScopeGuard>
#include <QTest>
#include <QTimer>

#include <iostream>

berry::IWorkbenchPage::Pointer QmitkDiffusionUiTest::GetPage()
{
  auto workbench = berry::PlatformUI::GetWorkbench();
  auto window = workbench->GetActiveWorkbenchWindow();
  // Offscreen, no window has the focus, so there may be no active window.
  if (window.IsNull() && !workbench->GetWorkbenchWindows().isEmpty())
    window = workbench->GetWorkbenchWindows().front();
  CPPUNIT_ASSERT_MESSAGE("No workbench window", window.IsNotNull());
  auto page = window->GetActivePage();
  CPPUNIT_ASSERT_MESSAGE("No workbench page", page.IsNotNull());
  return page;
}

QWidget *QmitkDiffusionUiTest::ShowView(const QString &viewId, const QString &controlsName)
{
  auto part = GetPage()->ShowView(viewId);
  CPPUNIT_ASSERT_MESSAGE("View " + viewId.toStdString() + " did not open", part.IsNotNull());
  Wait(50);

  for (auto *widget : QApplication::allWidgets())
  {
    if (widget->objectName() == controlsName)
      return widget;
  }
  CPPUNIT_FAIL("View " + viewId.toStdString() + " has no widget " + controlsName.toStdString());
  return nullptr;
}

void QmitkDiffusionUiTest::HideView(const QString &viewId)
{
  auto page = GetPage();
  auto part = page->FindView(viewId);
  if (part.IsNotNull())
    page->HideView(part);
  Wait(50);
}

bool QmitkDiffusionUiTest::IsViewOpen(const QString &viewId)
{
  return GetPage()->FindView(viewId).IsNotNull();
}

mitk::DataStorage::Pointer QmitkDiffusionUiTest::GetDataStorage()
{
  mitk::CoreServicePointer<mitk::IDataStorageService> service(mitk::CoreServices::GetDataStorageService());
  auto storage = service->GetActiveDataStorage();
  CPPUNIT_ASSERT_MESSAGE("No data storage", storage.IsNotNull());
  return storage;
}

void QmitkDiffusionUiTest::ClearDataStorage()
{
  auto storage = GetDataStorage();
  storage->Remove(storage->GetAll());
  Wait(50);
}

mitk::DataNode::Pointer QmitkDiffusionUiTest::Load(const std::string &testDataFile)
{
  auto nodes = mitk::IOUtil::Load(mitk::GetDiffusionTestDataFilePath(testDataFile), *GetDataStorage());
  CPPUNIT_ASSERT_MESSAGE("Loading " + testDataFile, nodes.IsNotNull() && nodes->Size() == 1);
  Wait(50);
  return nodes->front();
}

mitk::DataNode::Pointer QmitkDiffusionUiTest::AddNode(mitk::BaseData *data, const std::string &name, mitk::DataNode *parent)
{
  auto node = mitk::DataNode::New();
  node->SetData(data);
  node->SetName(name);
  GetDataStorage()->Add(node, parent);
  Wait(50);
  return node;
}

std::vector<mitk::DataNode::Pointer> QmitkDiffusionUiTest::GetNodesOfClass(const std::string &className)
{
  std::vector<mitk::DataNode::Pointer> result;
  const auto nodes = GetDataStorage()->GetAll();
  for (const auto &node : *nodes)
  {
    if (node->GetData() != nullptr && className == node->GetData()->GetNameOfClass())
      result.push_back(node);
  }
  return result;
}

void QmitkDiffusionUiTest::ChooseNode(QWidget *root, const char *boxName, mitk::DataNode *node)
{
  auto *box = GetChild<QmitkDataStorageComboBox>(root, boxName);
  // The box with a "none" entry has a SetSelectedNode of its own, which the base class does not
  // declare virtual.
  if (auto *boxWithNone = qobject_cast<QmitkDataStorageComboBoxWithSelectNone *>(box))
    boxWithNone->SetSelectedNode(node);
  else
    box->SetSelectedNode(node);
  Wait(20);
  CPPUNIT_ASSERT_MESSAGE(std::string(boxName) + " does not hold " + (node != nullptr ? node->GetName() : "none"),
                         box->GetSelectedNode().GetPointer() == node);
}

void QmitkDiffusionUiTest::Select(const QList<mitk::DataNode::Pointer> &nodes)
{
  auto page = GetPage();
  auto previous = page->GetActivePart();

  auto part = page->ShowView(QmitkDiffusionUiTestSelectionView::VIEW_ID);
  auto *view = dynamic_cast<QmitkDiffusionUiTestSelectionView *>(part.GetPointer());
  CPPUNIT_ASSERT_MESSAGE("No selection view", view != nullptr);
  view->Select(nodes);
  Wait(20);

  // The selection view shares a folder with the tested view; bringing the tested view back to the
  // front keeps its widgets visible, as a user sees them next to the Data Manager. The selection
  // stays: a view ignores its own (empty) selection.
  if (previous.IsNotNull() && previous != part)
    page->Activate(previous);
  Wait(20);
}

void QmitkDiffusionUiTest::InitializeViews()
{
  mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(GetDataStorage());
  Wait(50);
}

void QmitkDiffusionUiTest::SetCrosshair(const mitk::Point3D &position)
{
  auto *renderWindowPart = mitk::WorkbenchUtil::GetRenderWindowPart(GetPage(), mitk::WorkbenchUtil::OPEN);
  CPPUNIT_ASSERT_MESSAGE("No render window part", renderWindowPart != nullptr);
  renderWindowPart->SetSelectedPosition(position);
  Wait(100);
}

bool QmitkDiffusionUiTest::WaitFor(const std::function<bool()> &condition, int timeoutMs)
{
  return QTest::qWaitFor(condition, timeoutMs);
}

void QmitkDiffusionUiTest::Wait(int ms)
{
  QTest::qWait(ms);
}

std::string QmitkDiffusionUiTest::CreateTempDir()
{
  return mitk::IOUtil::CreateTemporaryDirectory("DiffusionUiTest_XXXXXX");
}

QmitkDiffusionUiTest::DialogHandler &QmitkDiffusionUiTest::DialogHandler::Instance()
{
  static DialogHandler handler;
  return handler;
}

void QmitkDiffusionUiTest::DialogHandler::Start()
{
  if (m_Timer != nullptr)
    return;
  m_Timer = new QTimer;
  m_Timer->setInterval(25);
  QObject::connect(m_Timer, &QTimer::timeout, m_Timer, [this]() { this->HandleModalWidget(); });
  m_Timer->start();
}

void QmitkDiffusionUiTest::DialogHandler::Stop()
{
  delete m_Timer;
  m_Timer = nullptr;
  m_Answers.clear();
}

void QmitkDiffusionUiTest::DialogHandler::AnswerMessageBox(QMessageBox::StandardButton button)
{
  m_Answers.push_back([button](QWidget *widget) {
    auto *box = qobject_cast<QMessageBox *>(widget);
    if (box == nullptr)
      return false;
    if (auto *buttonWidget = box->button(button))
      buttonWidget->click();
    else
      box->done(button);
    return true;
  });
}

void QmitkDiffusionUiTest::DialogHandler::AnswerFileDialog(const QString &path)
{
  m_Answers.push_back([path](QWidget *widget) {
    auto *dialog = qobject_cast<QFileDialog *>(widget);
    if (dialog == nullptr)
      return false;

    // The dialog's file system model loads folders asynchronously and may replace a selection
    // made before, and selectFile() leaves the name box alone while it has the focus; the answer
    // enters the folder, types the name and waits until the dialog reports the path. A folder is
    // chosen with an empty name, which the dialog takes as its current folder.
    const QFileInfo info(path);
    const QStringList expected = {info.absoluteFilePath()};
    const bool isFolder = dialog->fileMode() == QFileDialog::Directory && info.isDir();
    auto *nameEdit = dialog->findChild<QLineEdit *>("fileNameEdit");
    auto select = [dialog, nameEdit, &info, isFolder]() {
      dialog->setDirectory(isFolder ? info.absoluteFilePath() : info.absolutePath());
      if (nameEdit != nullptr)
        nameEdit->setText(isFolder ? QString() : info.fileName());
    };
    select();
    if (!QTest::qWaitFor([&]() { select(); return dialog->selectedFiles() == expected; }, 5000))
    {
      MITK_ERROR << "Diffusion UI test could not choose " << path.toStdString() << " in the file dialog (it has "
                 << dialog->selectedFiles().join(", ").toStdString() << ")";
      dialog->reject();
      return true;
    }
    // QFileDialog::accept() is protected; it is QDialog's public slot.
    static_cast<QDialog *>(dialog)->accept();
    return true;
  });
}

void QmitkDiffusionUiTest::DialogHandler::AnswerColorDialog(const QColor &color)
{
  m_Answers.push_back([color](QWidget *widget) {
    auto *dialog = qobject_cast<QColorDialog *>(widget);
    if (dialog == nullptr)
      return false;
    dialog->setCurrentColor(color);
    dialog->accept();
    return true;
  });
}

QStringList QmitkDiffusionUiTest::DialogHandler::TakeMessages()
{
  QStringList messages;
  messages.swap(m_Messages);
  return messages;
}

QStringList QmitkDiffusionUiTest::DialogHandler::PeekMessages() const
{
  return m_Messages;
}

void QmitkDiffusionUiTest::DialogHandler::HandleModalWidget()
{
  // An answer may process events (file dialog), which runs this timer again.
  if (m_Busy)
    return;
  auto *widget = QApplication::activeModalWidget();
  if (widget == nullptr)
    return;
  m_Busy = true;
  const auto done = qScopeGuard([this]() { m_Busy = false; });

  auto *box = qobject_cast<QMessageBox *>(widget);
  if (box != nullptr)
  {
    m_Messages.push_back(box->windowTitle() + ": " + box->text());
    MITK_INFO << "Diffusion UI test closes message box \"" << box->windowTitle().toStdString() << "\": "
              << box->text().toStdString();
  }

  // An answer is for the dialog type it was queued for; another dialog in between is cancelled.
  if (!m_Answers.empty() && m_Answers.front()(widget))
  {
    m_Answers.pop_front();
    return;
  }

  if (box != nullptr)
  {
    if (auto *escape = box->escapeButton())
      escape->click();
    else
      box->reject();
  }
  else if (auto *dialog = qobject_cast<QDialog *>(widget))
  {
    MITK_INFO << "Diffusion UI test cancels dialog " << widget->metaObject()->className();
    dialog->reject();
  }
  else
  {
    widget->close();
  }
}

QmitkDiffusionUiTest::LogCapture::LogCapture()
{
  mitk::RegisterBackend(this);
}

QmitkDiffusionUiTest::LogCapture::~LogCapture()
{
  mitk::UnregisterBackend(this);
}

void QmitkDiffusionUiTest::LogCapture::ProcessMessage(const mitk::LogMessage &message)
{
  std::lock_guard<std::mutex> lock(m_Mutex);
  m_Messages.push_back(message.Message);
  if (message.Level == mitk::LogLevel::Error || message.Level == mitk::LogLevel::Fatal)
  {
    m_Errors.push_back(message.Message);
    std::cout << "Diffusion UI test: logged error: " << message.Message << std::endl;
  }
}

mitk::LogBackendBase::OutputType QmitkDiffusionUiTest::LogCapture::GetOutputType() const
{
  return OutputType::Other;
}

bool QmitkDiffusionUiTest::LogCapture::Contains(const std::string &text) const
{
  std::lock_guard<std::mutex> lock(m_Mutex);
  for (const auto &message : m_Messages)
  {
    if (message.find(text) != std::string::npos)
      return true;
  }
  return false;
}

bool QmitkDiffusionUiTest::LogCapture::ContainsError() const
{
  std::lock_guard<std::mutex> lock(m_Mutex);
  return !m_Errors.empty();
}

std::vector<std::string> QmitkDiffusionUiTest::LogCapture::GetMessages() const
{
  std::lock_guard<std::mutex> lock(m_Mutex);
  return m_Messages;
}
