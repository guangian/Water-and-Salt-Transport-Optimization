import { Navigate, Route, Routes } from 'react-router-dom';
import { AppShell } from './components/layout/AppShell';
import DashboardPage from './pages/Dashboard';
import NewProjectPage from './pages/NewProject';
import ProjectPage from './pages/Project';
import SettingsPage from './pages/Settings';
import KnowledgePage from './pages/Knowledge';

export default function App() {
  return (
    <AppShell>
      <Routes>
        <Route path="/" element={<DashboardPage />} />
        <Route path="/new" element={<NewProjectPage />} />
        <Route path="/projects/:id" element={<ProjectPage />} />
        <Route path="/settings" element={<SettingsPage />} />
        <Route path="/knowledge" element={<KnowledgePage />} />
        <Route path="/knowledge/:slug" element={<KnowledgePage />} />
        <Route path="*" element={<Navigate to="/" replace />} />
      </Routes>
    </AppShell>
  );
}
