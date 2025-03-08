Sub SummaryOfConsMarkers()
    Dim ws As Worksheet
    Dim SummarySheet As Worksheet
    Dim SummaryWB As Workbook
    Dim UniqueStrings As Collection
    Dim UniqueString As Variant
    Dim CellValue As Variant
    Dim Header As Range
    Dim StringCount As Long
    Dim LastRow As Long
    Dim Headers() As String
    Dim SheetIndex As Integer
    Dim ColIndex As Integer
    Dim RowIndex As Integer

    Application.ScreenUpdating = False
    Application.Calculation = xlCalculationManual
    
    Set SummaryWB = Workbooks.Add
    Set SummarySheet = SummaryWB.Sheets(1)
    SummarySheet.Name = "Summary"
    
    ' Create headers for the summary table with all unique strings
    Set UniqueStrings = New Collection
    On Error Resume Next
    For Each ws In ThisWorkbook.Sheets
        LastRow = ws.Cells(ws.Rows.Count, "N").End(xlUp).Row
        For Each CellValue In ws.Range("N1:N" & LastRow)
            ' Add only non-empty and non-numeric strings as unique values
            If Not IsEmpty(CellValue) And Not IsNumeric(CellValue) Then
                UniqueStrings.Add CellValue.Value, CStr(CellValue.Value)
            End If
        Next CellValue
    Next ws
    On Error GoTo 0 ' Reset error handling
    
    ' Create header row
    ColIndex = 2
    For Each UniqueString In UniqueStrings
        SummarySheet.Cells(1, ColIndex).Value = UniqueString
        ColIndex = ColIndex + 1
    Next UniqueString
    ReDim Headers(1 To ColIndex - 2)
    For i = 1 To UBound(Headers)
        Headers(i) = SummarySheet.Cells(1, i + 1).Value
    Next i
    
    ' Fill in the counts for each unique string
    SheetIndex = 2
    For Each ws In ThisWorkbook.Sheets
        SummarySheet.Cells(SheetIndex, 1).Value = ws.Name
        For i = 1 To UBound(Headers)
            StringCount = Application.WorksheetFunction.CountIf(ws.Range("N:N"), Headers(i))
            SummarySheet.Cells(SheetIndex, i + 1).Value = StringCount
        Next i
        SheetIndex = SheetIndex + 1
    Next ws
    
    SummarySheet.Columns.AutoFit
    
    Application.ScreenUpdating = True
    Application.Calculation = xlCalculationAutomatic
    
    MsgBox "Summary table has been created.", vbInformation
End Sub
