Sub SearchForMarkers()
    Dim ws As Worksheet
    Dim SearchColumn As String
    Dim CheckColumn As String
    Dim OutputColumn As String
    Dim Cell As Range
    Dim FirstAddress As String
    Dim Category As Variant
    Dim Markers As Collection
    Dim EndoMarkers As Collection
    Dim EpiMarkers As Collection
    Dim MyoMarkers As Collection
    Dim MesMarkers As Collection
    Dim MPMarkers As Collection
    Dim HepMarkers As Collection
    Dim MarkerGroup As Collection
    Dim i As Integer

    ' Initialize marker collections
    Set EndoMarkers = New Collection
    Set EpiMarkers = New Collection
    Set MyoMarkers = New Collection
    Set MesMarkers = New Collection
    Set MPMarkers = New Collection
    Set HepMarkers = New Collection

    ' Add markers to each collection
    With EndoMarkers
        .Add "Emcn": .Add "Cdh5": .Add "Npr3": .Add "Pecam1": .Add "Egfl7"
        .Add "Plvap": .Add "Ecscr": .Add "Icam2": .Add "Klf2": .Add "Plxnd1"
    End With
    With EpiMarkers
        .Add "Upk3b": .Add "Tbx18": .Add "Wt1": .Add "Aldh1a2": .Add "Upk1b"
        .Add "Sparc": .Add "Kcne1l": .Add "Tmem255a"
    End With
    With MyoMarkers
        .Add "Actc1": .Add "Ttn": .Add "Myl4": .Add "Tnnc1": .Add "Tnnt2"
        .Add "Acta2": .Add "Nebl": .Add "Cryab"
    End With
    With MesMarkers
        .Add "Postn": .Add "Cthrc1": .Add "Sox9": .Add "Pdgfra": .Add "Papss2"
    End With
    With MPMarkers
        .Add "Osr1": .Add "Rgs5": .Add "Foxf1": .Add "Isl1": .Add "Tbx1"
        .Add "Fgf10": .Add "Mfap4"
    End With
    With HepMarkers
        .Add "Alb": .Add "Afp": .Add "Apoa1": .Add "Apom": .Add "Hnf4a"
    End With

    ' Create a collection of collections
    Set Markers = New Collection
    Markers.Add EndoMarkers: Markers.Add EpiMarkers
    Markers.Add MyoMarkers: Markers.Add MesMarkers
    Markers.Add MPMarkers: Markers.Add HepMarkers

    ' Define the columns to search and modify
    SearchColumn = "A" ' Column where you are searching for the string
    CheckColumn1 = "C"  ' First column where you are checking if the value is positive
    CheckColumn2 = "H"  ' Second column where you are checking if the value is positive
    OutputColumn = "N" ' Column where you are writing the category of the string

    ' Loop through each sheet in the workbook
    For Each ws In ThisWorkbook.Worksheets
        ' Loop through each collection of markers
        For i = 1 To Markers.Count
            Set MarkerGroup = Markers(i)
            ' Loop through each marker in the collection
            For Each Category In MarkerGroup
                ' Use Find method to search for the string with exact match
                Set Cell = ws.Columns(SearchColumn).Find(What:=Category, LookIn:=xlValues, LookAt:=xlWhole)
                FirstAddress = ""
                Do While Not Cell Is Nothing
                    ' Save the address of the first found cell
                    If FirstAddress = "" Then
                        FirstAddress = Cell.Address
                    ElseIf Cell.Address = FirstAddress Then
                        ' If we're back at the first cell, exit the loop
                        Exit Do
                    End If

                    ' Check if the corresponding cell in CheckColumn is positive
                    If ws.Cells(Cell.Row, CheckColumn1).Value > 0 And _
                       (ws.Cells(Cell.Row, CheckColumn2).Value > 0 Or IsEmpty(ws.Cells(Cell.Row, CheckColumn2).Value)) Then
                        ' Change format if positive or column H is blank
                        With Cell.Font
                            .Bold = True
                            .Color = RGB(255, 0, 0) ' Red color
                        End With
                        ' Write the category to OutputColumn based on the collection
                        Select Case i
                            Case 1: ws.Cells(Cell.Row, OutputColumn).Value = "Endocardial"
                            Case 2: ws.Cells(Cell.Row, OutputColumn).Value = "Epicardial"
                            Case 3: ws.Cells(Cell.Row, OutputColumn).Value = "Myocardial"
                            Case 4: ws.Cells(Cell.Row, OutputColumn).Value = "Mesenchymal"
                            Case 5: ws.Cells(Cell.Row, OutputColumn).Value = "MP"
                            Case 6: ws.Cells(Cell.Row, OutputColumn).Value = "Hepatocyte"
                        End Select
                    End If
                    
                    ' Find the next cell with the exact match
                    Set Cell = ws.Columns(SearchColumn).Find(What:=Category, After:=Cell, LookIn:=xlValues, LookAt:=xlWhole)
                Loop
                FirstAddress = ""
            Next Category
        Next i
    Next ws
End Sub

