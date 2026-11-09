/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef QmitkDiffusionWarningText_h
#define QmitkDiffusionWarningText_h

#include <QmitkIconTheme.h>

#include <QString>

/**
 * \brief Warning texts of the diffusion views in the colour of the current theme.
 *
 * The themes colour warnings with a "font.warning" rule, which reaches rich text only through a
 * QTextDocument given the application style sheet; a QLabel does not apply it to its text. The
 * views' warning labels therefore take the colour the rule declares (QmitkIconTheme) as an
 * explicit font colour. The class stays on the tag, so the text still reads as a warning where
 * the rule does apply.
 */
namespace QmitkDiffusionWarningText
{
  /** \brief Rich text showing text as a warning: bold, in the theme's warning colour if it has one. */
  inline QString Format(const QString &text)
  {
    const auto color = QmitkIconTheme::GetWarningTextColor();
    return color.isEmpty()
      ? QStringLiteral("<font class=\"warning\"><b>%1</b></font>").arg(text)
      : QStringLiteral("<font class=\"warning\" color=\"%1\"><b>%2</b></font>").arg(color, text);
  }
}

#endif
